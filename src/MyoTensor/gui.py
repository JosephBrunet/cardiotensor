# -*- coding: utf-8 -*-
import sys
import cv2
import numpy as np
import argparse
from pathlib import Path
from PyQt5.QtGui import (
    QBrush, QPainter, QPen, QPixmap, QKeySequence, QColor, QImage,
    QDoubleValidator, QIntValidator
)
from PyQt5.QtWidgets import (
    QFileDialog, QApplication, QGraphicsEllipseItem, QGraphicsItem,
    QGraphicsLineItem, QGraphicsScene, QGraphicsView, QHBoxLayout,
    QPushButton, QVBoxLayout, QWidget, QShortcut, QLineEdit, QLabel,
    QMenuBar, QAction, QSpinBox, QRadioButton, QButtonGroup
)
from skimage import color, io
from skimage.measure import block_reduce
from MyoTensor.utils import *
from MyoTensor.processing_functions import *
from MyoTensor.analysis_functions import *
import matplotlib.pyplot as plt



def np2pixmap(np_img):
    height, width, channel = np_img.shape
    bytes_per_line = 3 * width
    q_img = QImage(np_img.data, width, height, bytes_per_line, QImage.Format_RGB888)
    return QPixmap.fromImage(q_img)


class Window(QWidget):
    def __init__(self, conf_file_path, N_slice, N_line, angle_range, image_mode):
        super().__init__()

        self.half_point_size = 5
        self.line_np = None
        self.color_idx = 0
        self.bg_img = None
        self.is_mouse_down = False
        self.point_size = self.half_point_size * 2
        self.start_point = None
        self.end_point = None
        self.start_pos = (None, None)
        
        self.N_slice = N_slice
        self.angle_range = angle_range
        self.N_line = N_line
        self.image_mode = image_mode
        self.intensity_profiles = []
        self.x_min_lim = None
        self.x_max_lim = None
        self.y_min_lim = None
        self.y_max_lim = None

        self.view = QGraphicsView()
        self.view.setRenderHint(QPainter.Antialiasing)


        params = read_conf_file(conf_file_path)
        try:
            params = read_conf_file(conf_file_path)
        except Exception as e:
            sys.exit(f'⚠️  Error reading parameter file: {conf_file_path}')

        VOLUME_PATH, MASK_PATH, IS_FLIP, OUTPUT_DIR, OUTPUT_TYPE, SIGMA, RHO, N_CHUNK, PT_MV, PT_APEX, IS_TEST, N_SLICE_TEST = [
            params[key] for key in [
                'IMAGES_PATH', 'MASK_PATH', 'FLIP', 'OUTPUT_PATH', 'OUTPUT_TYPE',
                'SIGMA', 'RHO', 'N_CHUNK', 'POINT_MITRAL_VALVE', 'POINT_APEX', 'TEST', 'N_SLICE_TEST'
            ]
        ]
        
        self.MASK_PATH = MASK_PATH
        self.OUTPUT_DIR = OUTPUT_DIR
        
        HA_path = Path(self.OUTPUT_DIR) / "HA"
        IA_path = Path(self.OUTPUT_DIR) / "IA"
        FA_path = Path(self.OUTPUT_DIR) / "FA"
        
        if not HA_path.exists:
            sys.exit(f"No HA folder ({HA_path})")
        if not IA_path.exists:
            sys.exit(f"No IA folder ({IA_path})")
        elif not FA_path.exists:
            sys.exit(f"No FA folder ({FA_path})")


        # Select the appropriate image list based on the output mode
        if self.image_mode == 'HA':
            self.img_list, img_type = get_image_list(HA_path)
        elif self.image_mode == 'IA':
            self.img_list, img_type = get_image_list(IA_path)
        elif self.image_mode == 'FA':
            self.img_list, img_type = get_image_list(FA_path)
        else:
            sys.exit(f"Invalid output mode: {self.image_mode}. Must be 'HA', 'IA', or 'FA'.")

        # Compute center line and center point
        center_line = interpolate_points(PT_MV, PT_APEX, len(self.img_list))
        center_point = np.around(center_line[self.N_slice])

        # Load the selected slice
        slice_path = self.img_list[self.N_slice]
        if not slice_path.exists():
            print(f"File {slice_path} doesn't exist")
            return
        else:
            print(f"Slice found ({slice_path})")
            
            
        self.current_img = self.load_image(slice_path, self.MASK_PATH)


        self.load_image_to_gui(self.current_img)
        
        menu_bar = QMenuBar(self)
        file_menu = menu_bar.addMenu('File')
        
        change_slice_action = QAction('Change Slice', self)
        change_slice_action.triggered.connect(self.change_slice_dialog)
        file_menu.addAction(change_slice_action)
        
        change_mode = QAction('Choose image type', self)
        change_mode.triggered.connect(self.change_img_mode)
        file_menu.addAction(change_mode)

        graph_menu = menu_bar.addMenu('Graph')
        set_x_lim_action = QAction('Set x axis limits', self)
        set_x_lim_action.triggered.connect(self.set_x_lim_dialog)
        graph_menu.addAction(set_x_lim_action)
        
        set_y_lim_action = QAction('Set y axis limits', self)
        set_y_lim_action.triggered.connect(self.set_y_lim_dialog)
        graph_menu.addAction(set_y_lim_action)

        vbox = QVBoxLayout(self)
        vbox.setMenuBar(menu_bar)
        vbox.addWidget(self.view)

        plot_profile_button = QPushButton("Plot profile")
        save_profile_button = QPushButton("Save Profile")

        label_angle_range = QLabel("Angle range:")
        self.input_angle_range = QLineEdit(self)
        self.input_angle_range.setValidator(QDoubleValidator(0, 360, 2))
        self.input_angle_range.setText(str(self.angle_range))
        
        label_N_line = QLabel("Number of lines:")
        self.input_N_line = QLineEdit(self)
        self.input_N_line.setValidator(QIntValidator(1, 9999))
        self.input_N_line.setText(str(self.N_line))

        input_layout = QHBoxLayout()
        input_layout.addWidget(label_angle_range)
        input_layout.addWidget(self.input_angle_range)
        input_layout.addWidget(label_N_line)
        input_layout.addWidget(self.input_N_line)

        hbox = QHBoxLayout(self)
        hbox.addWidget(plot_profile_button)
        hbox.addWidget(save_profile_button)
        
        vbox.addLayout(input_layout)
        vbox.addLayout(hbox)

        self.setLayout(vbox)
        
        self.input_angle_range.textChanged.connect(self.update_text)
        self.input_N_line.textChanged.connect(self.update_text)

        self.quit_shortcut = QShortcut(QKeySequence("Ctrl+Q"), self)
        self.quit_shortcut.activated.connect(lambda: quit())

        # self.undo_shortcut = QShortcut(QKeySequence("Ctrl+Z"), self)
        # self.undo_shortcut.activated.connect(self.undo)
        
        plot_profile_button.clicked.connect(self.plot_profile)
        save_profile_button.clicked.connect(self.save_profile)

    # def undo(self):
    #     if self.prev_mask is None:
    #         print("No previous mask record")
    #         return

    #     self.color_idx -= 1

    #     bg = Image.fromarray(self.current_img_rgb.astype("uint8"), "RGB")
    #     mask = Image.fromarray(self.prev_mask.astype("uint8"), "RGB")
    #     img = Image.blend(bg, mask, 0.2)

    #     self.scene.removeItem(self.bg_img)
    #     self.bg_img = self.scene.addPixmap(np2pixmap(np.array(img)))

    #     self.mask_c = self.prev_mask
    #     self.prev_mask = None

    def update_text(self):
        angle_range_text = self.input_angle_range.text()
        if not angle_range_text or float(angle_range_text) <= 0:
            self.angle_range = 0.1
            self.input_angle_range.setText("0.1")
        else:
            self.angle_range = float(angle_range_text)
                
        n_line_text = self.input_N_line.text()
        if not n_line_text or int(n_line_text) < 1:
            self.N_line = 1
            self.input_N_line.setText("1")
        else:
            self.N_line = int(n_line_text)
            
        self.update_plot()
        
    def update_plot(self):
        if not self.start_point or not self.end_point:
            return

        # Clear existing lines from the scene
        
        if self.lines:
            for l in self.lines:
                if l.scene() == self.scene:
                    self.scene.removeItem(l)
        self.lines = []

        # Recalculate end points based on the updated angle range and number of lines
        sx, sy = self.start_pos
        ex, ey = self.end_point.rect().center().x(), self.end_point.rect().center().y()

        self.end_points = find_end_points([sy, sx], [ey, ex], float(self.angle_range), int(self.N_line))

        # Draw new lines
        for i in range(self.N_line):
            line = self.scene.addLine(sx, sy, self.end_points[i][1], self.end_points[i][0], pen=QPen(QColor("black"), 2))
            self.lines.append(line)

        # Redraw the start and end points
        if self.end_point is not None and self.end_point.scene() == self.scene:
            self.scene.removeItem(self.end_point)
        if self.start_point is not None and self.start_point.scene() == self.scene:
            self.scene.removeItem(self.start_point)

        self.start_point = self.scene.addEllipse(
            sx - self.half_point_size,
            sy - self.half_point_size,
            self.point_size,
            self.point_size,
            pen=QPen(QColor("black")),
            brush=QBrush(QColor("black"))
        )

        self.end_point = self.scene.addEllipse(
            ex - self.half_point_size,
            ey - self.half_point_size,
            self.point_size,
            self.point_size,
            pen=QPen(QColor("black")),
            brush=QBrush(QColor("black"))
        )
        

    def plot_profile(self):
        if self.line_np is None:
            print("No line drawn")
            return
        if not self.line_np.any():
            print("No line drawn")
            return
        
        start_point = self.line_np[0:2] * self.bin_factor
        end_point = self.line_np[2:] * self.bin_factor
        
        print(start_point)
        print(end_point)
                
        # Adjust minimum and maximum based on the mode
        if self.image_mode == 'HA':
            minimum = -90  # Use the default behavior or adjust values
            maximum = 90
            label_y = 'Helical Angle (°)'
        elif self.image_mode == 'IA':
            # Adjust min/max for IA if necessary
            minimum = -90
            maximum = 90
            label_y = 'Intrusion Angle (°)'
        elif self.image_mode == 'FA':
            # Adjust min/max for FA if necessary
            minimum = 0
            maximum = 1
            label_y = 'Fractional Anisotropy'
                
        self.intensity_profiles = calculate_intensities(self.current_img, start_point, end_point, self.angle_range, self.N_line,  max_value = maximum, min_value = minimum)
        
        plot_intensity(self.intensity_profiles, label_y=label_y, x_max_lim=self.x_max_lim, x_min_lim=self.x_min_lim, y_max_lim=self.y_max_lim, y_min_lim=self.y_min_lim)

        

        # bin_factor = 16                      
        # current_img_bin = block_reduce(self.current_img, block_size=(bin_factor, bin_factor), func=np.mean)
        
        # intensity_profiles = calculate_intensities(current_img_bin, start_point/bin_factor, end_point/bin_factor, self.angle_range, self.N_line)      
        # plot_intensity(intensity_profiles)
        
        
        


    def save_profile(self):
        if self.line_np is None or not self.line_np.any():
            print("No line drawn")
            return
        
        start_point = self.line_np[0:2] * self.bin_factor
        end_point = self.line_np[2:] * self.bin_factor

        if not self.intensity_profiles:
            self.intensity_profiles = calculate_intensities(self.current_img, start_point, end_point, self.angle_range, self.N_line)
        
        save_path, _ = QFileDialog.getSaveFileName(self, "Save Profile", "", "Csv Files (*.csv);;All Files (*)")
        
        if save_path:
            if not save_path.lower().endswith('.csv'):
                save_path += '.csv'
            save_intensity(self.intensity_profiles, save_path)


    def load_image(self, file_path, mask_path='', bin_factor=None):
        img = cv2.imread(str(file_path), -1)
        if bin_factor:
            img = block_reduce(img, block_size=(bin_factor, bin_factor), func=np.mean)
        
        img64 = img.astype(float)
        
        if mask_path:
            print(f"\n---------------------------------")
            print(f"READING MASK INFORMATION FROM {mask_path}...\n")
            mask_list, mask_type = get_image_list(mask_path)
            N_mask = len(mask_list)
            mask_bin_factor = len(self.img_list) / N_mask
            print(f"Mask bining factor: {mask_bin_factor}\n")
            
            img_mask = self.load_image(mask_list[int(self.N_slice / mask_bin_factor)])
            img_mask = cv2.resize(img_mask, (img64.shape[1], img64.shape[0]), interpolation = cv2.INTER_LINEAR)
            
            assert img_mask.shape == img.shape, f"Mask shape {mask.shape} does not match volume shape {volume.shape}"           

            img64[img_mask == 0] = np.nan       
        
        return img64

        

    def load_image_to_gui(self, current_img):
                
        height, width = current_img.shape[:2]
        
        # Calculate the bin_factor such that the downsampled image is less than 1000 pixels in both dimensions
        self.bin_factor = 1
        while height // self.bin_factor >= 1000 or width // self.bin_factor >= 1000:
            self.bin_factor *= 2
        
        print(f"Bin factor: {self.bin_factor}")
                
        current_img_bin = block_reduce(current_img, block_size=(self.bin_factor, self.bin_factor), func=np.mean)
        
        
        # current_img_bin = current_img_bin.astype(float)

        
        minimum = np.nanmin(current_img_bin)
        maximum = np.nanmax(current_img_bin)
        current_img_bin = (current_img_bin + np.abs(minimum)) * (1 / (maximum - minimum)) 
        cmap = plt.get_cmap('hsv')
        current_img_rgb = cmap(current_img_bin)
        current_img_rgb = (current_img_rgb[:, :, :3] * 255).astype(np.uint8)

        gray_color = [128, 128, 128]
        current_img_rgb[np.isnan(current_img_bin)] = gray_color
        
        self.current_img_rgb = current_img_rgb

        pixmap = np2pixmap(self.current_img_rgb)

        H, W, _ = self.current_img_rgb.shape

        self.scene = QGraphicsScene(0, 0, W, H)
        self.end_point = None
        self.start_point = None
        self.line = None
        self.lines = []
        self.bg_img = self.scene.addPixmap(pixmap)
        self.bg_img.setPos(0, 0)
        
        # Set the scene rectangle to the size of the image
        self.scene.setSceneRect(0, 0, W, H)
        self.view.setScene(self.scene)
        
        
        self.scene.mousePressEvent = self.mouse_press
        self.scene.mouseMoveEvent = self.mouse_move
        self.scene.mouseReleaseEvent = self.mouse_release


    def mouse_press(self, ev):
        x, y = ev.scenePos().x(), ev.scenePos().y()
        self.is_mouse_down = True
        self.start_pos = ev.scenePos().x(), ev.scenePos().y()
        
        try:
            if self.start_point is not None:
                self.scene.removeItem(self.start_point)
        except:
            pass
                
        self.start_point = self.scene.addEllipse(
            x - self.half_point_size,
            y - self.half_point_size,
            self.point_size,
            self.point_size,
            pen=QPen(QColor("black")),
            brush=QBrush(QColor("black")),
        )

    def mouse_move(self, ev):
        if not self.is_mouse_down:
            return

        x, y = ev.scenePos().x(), ev.scenePos().y()

        if self.end_point is not None and self.end_point.scene() == self.scene:
            self.scene.removeItem(self.end_point)
            
        self.end_point = self.scene.addEllipse(
            x - self.half_point_size,
            y - self.half_point_size,
            self.point_size,
            self.point_size,
            pen=QPen(QColor("black")),
            brush=QBrush(QColor("black")),
        )

        # if self.line is not None:
        #     self.scene.removeItem(self.line)
        # self.line = self.scene.addLine(
        #     sx, sy, x, y, pen=QPen(QColor("black"))
        # )
        
        
        sx, sy = self.start_pos
        self.end_points = find_end_points([sy, sx], [y, x], float(self.angle_range), int(self.N_line))

        if self.lines:
            for l in self.lines:
                if l.scene() == self.scene:
                    self.scene.removeItem(l)
                
        for i in range(self.N_line):
            self.lines.append(self.scene.addLine(
                sx, sy, self.end_points[i][1], self.end_points[i][0], pen=QPen(QColor("black"), 2)
            ))

    def mouse_release(self, ev):
        x, y = ev.scenePos().x(), ev.scenePos().y()
        sx, sy = self.start_pos

        self.is_mouse_down = False

        H, W, _ = self.current_img_rgb.shape
        self.line_np = np.array([sy, sx, y, x])
        print("Line:", self.line_np * self.bin_factor)


    def change_slice_dialog(self):
        self.dialog = QWidget()
        self.dialog.setWindowTitle("Change Slice")

        layout = QVBoxLayout()

        self.spin_box = QSpinBox()
        self.spin_box.setMinimum(0)
        self.spin_box.setMaximum(len(self.img_list) - 1)
        self.spin_box.setValue(self.N_slice)
        layout.addWidget(self.spin_box)

        change_button = QPushButton("Change")
        change_button.clicked.connect(self.change_slice)
        layout.addWidget(change_button)

        self.dialog.setLayout(layout)
        self.dialog.show()
        
    def change_slice(self):
        self.N_slice = self.spin_box.value()
        slice_path = self.img_list[self.N_slice]
        self.current_img = self.load_image(slice_path, self.MASK_PATH).astype(float)
        self.load_image_to_gui(self.current_img)
        self.dialog.close()
        
        
        
                
    def set_x_lim_dialog(self):
        self.dialog = QWidget()
        self.dialog.setWindowTitle("Change X-Axis Limits")

        layout = QVBoxLayout()

        # Create spin boxes for min and max limits
        self.x_min_spin_box = QSpinBox()
        self.x_max_spin_box = QSpinBox()

        # Set ranges for spin boxes (can be modified as needed)
        self.x_min_spin_box.setRange(-10000, 10000)  # Set appropriate limits
        self.x_max_spin_box.setRange(-10000, 10000)

        # Add spin boxes to the layout
        layout.addWidget(QLabel("X Min:"))
        layout.addWidget(self.x_min_spin_box)
        layout.addWidget(QLabel("X Max:"))
        layout.addWidget(self.x_max_spin_box)

        # Create the change button
        change_button = QPushButton("Set X Limits")
        change_button.clicked.connect(self.set_x_lim)
        layout.addWidget(change_button)

        # Set layout and display dialog
        self.dialog.setLayout(layout)
        self.dialog.show()

    def set_x_lim(self):
        # Get values from spin boxes
        self.x_min_lim = self.x_min_spin_box.value()
        self.x_max_lim = self.x_max_spin_box.value()

        # Close the dialog after setting limits
        self.dialog.close()

    def set_y_lim_dialog(self):
        self.dialog = QWidget()
        self.dialog.setWindowTitle("Change Y-Axis Limits")

        layout = QVBoxLayout()

        # Create spin boxes for min and max limits
        self.y_min_spin_box = QSpinBox()
        self.y_max_spin_box = QSpinBox()

        # Set ranges for spin boxes (can be modified as needed)
        self.y_min_spin_box.setRange(-10000, 10000)  # Set appropriate limits
        self.y_max_spin_box.setRange(-10000, 10000)

        # Add spin boxes to the layout
        layout.addWidget(QLabel("Y Min:"))
        layout.addWidget(self.y_min_spin_box)
        layout.addWidget(QLabel("Y Max:"))
        layout.addWidget(self.y_max_spin_box)

        # Create the change button
        change_button = QPushButton("Set Y Limits")
        change_button.clicked.connect(self.set_y_lim)
        layout.addWidget(change_button)

        # Set layout and display dialog
        self.dialog.setLayout(layout)
        self.dialog.show()

    def set_y_lim(self):
        # Get values from spin boxes
        self.y_min_lim = self.y_min_spin_box.value()
        self.y_max_lim = self.y_max_spin_box.value()


        # Close the dialog after setting limits
        self.dialog.close()

            
            
         
            
    def change_img_mode(self):
        # Create a new dialog for changing the image mode
        self.dialog = QWidget()
        self.dialog.setWindowTitle("Choose Image Mode")

        layout = QVBoxLayout()

        # Create radio buttons for each mode (HA, IA, FA)
        ha_radio = QRadioButton("HA", self.dialog)
        ia_radio = QRadioButton("IA", self.dialog)
        fa_radio = QRadioButton("FA", self.dialog)

        # Create a button group to ensure only one is selected at a time
        button_group = QButtonGroup(self.dialog)
        button_group.addButton(ha_radio)
        button_group.addButton(ia_radio)
        button_group.addButton(fa_radio)

        # Set the currently active mode as checked
        if self.image_mode == "HA":
            ha_radio.setChecked(True)
        elif self.image_mode == "IA":
            ia_radio.setChecked(True)
        elif self.image_mode == "FA":
            fa_radio.setChecked(True)

        # Add radio buttons to the layout
        layout.addWidget(ha_radio)
        layout.addWidget(ia_radio)
        layout.addWidget(fa_radio)

        # Create a button to confirm the selection
        confirm_button = QPushButton("Confirm", self.dialog)
        confirm_button.clicked.connect(lambda: self.set_img_mode(ha_radio, ia_radio, fa_radio))

        # Add the confirm button to the layout
        layout.addWidget(confirm_button)

        # Set the layout for the dialog and display it
        self.dialog.setLayout(layout)
        self.dialog.show()

    def set_img_mode(self, ha_radio, ia_radio, fa_radio):
        # Set the image mode based on the selected radio button
        if ha_radio.isChecked():
            self.image_mode = "HA"
        elif ia_radio.isChecked():
            self.image_mode = "IA"
        elif fa_radio.isChecked():
            self.image_mode = "FA"

        print(f"Image mode changed to {self.image_mode}")

        # Reload the images based on the selected mode
        if self.image_mode == 'HA':
            self.img_list, img_type = get_image_list(Path(self.OUTPUT_DIR) / "HA")
        elif self.image_mode == 'IA':
            self.img_list, img_type = get_image_list(Path(self.OUTPUT_DIR) / "IA")
        elif self.image_mode == 'FA':
            self.img_list, img_type = get_image_list(Path(self.OUTPUT_DIR) / "FA")
        else:
            print("Invalid image mode")

        # Load the first slice in the new mode
        self.current_img = self.load_image(self.img_list[self.N_slice], self.MASK_PATH)
        self.load_image_to_gui(self.current_img)

        # Close the dialog
        self.dialog.close()
                
            
            
            
            
            
            
            
            

def parse_arguments():
    parser = argparse.ArgumentParser(description="Image Processing Script")
    parser.add_argument("conf_file_path", type=str, help="Path to the configuration file")
    parser.add_argument("N_slice", type=int, help="Slice number")
    parser.add_argument("--N_line", type=int, default=5, help="Number of lines")
    parser.add_argument("--angle_range", type=float, default=20, help="Angle range in degrees")
    parser.add_argument("--image_mode", type=str, default="HA", help="Output mode (HA, IA, or FA)")

    return parser.parse_args()

def main():
    args = parse_arguments()

    app = QApplication(sys.argv)
    w = Window(args.conf_file_path, args.N_slice, args.N_line, args.angle_range, args.image_mode)
    w.show()
    app.exec()


if __name__ == "__main__":
    main()