import sys
from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Colormap
from PyQt5.QtCore import Qt
from PyQt5.QtGui import (
    QBrush,
    QColor,
    QDoubleValidator,
    QImage,
    QIntValidator,
    QKeySequence,
    QPainter,
    QPen,
    QPixmap,
)
from PyQt5.QtWidgets import (
    QAction,
    QButtonGroup,
    QFileDialog,
    QGraphicsEllipseItem,
    QGraphicsLineItem,
    QGraphicsPixmapItem,
    QGraphicsScene,
    QGraphicsSceneMouseEvent,
    QGraphicsView,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenuBar,
    QPushButton,
    QRadioButton,
    QShortcut,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)
from skimage.measure import block_reduce

from cardiotensor.analysis.analysis_functions import (
    calculate_intensities,
    find_end_points,
    plot_intensity,
    save_intensity,
)
from cardiotensor.colormaps.helix_angle import helix_angle_cmap
from cardiotensor.utils.DataReader import DataReader

# ---------- small helpers ----------


def np2pixmap(np_img: np.ndarray) -> QPixmap:
    h, w, _ = np_img.shape
    q_img = QImage(np_img.data, w, h, 3 * w, QImage.Format_RGB888)
    return QPixmap.fromImage(q_img)


def discover_modes(base: Path) -> list[str]:
    """Return list of available subfolders among HA, IA, AZ, EL, FA, in a stable order."""
    order = ["HA", "IA", "EL", "AZ", "FA"]
    return [m for m in order if (base / m).exists() and (base / m).is_dir()]


def default_cmap_for(mode: str):
    m = mode.upper()
    if m in {"HA", "IA", "EL"}:
        return helix_angle_cmap
    if m == "AZ":
        return plt.get_cmap("hsv")
    return plt.get_cmap("inferno")  # FA or fallback


def convert_slice_for_display(slice2d: np.ndarray, mode: str) -> np.ndarray:
    """Return normalized float slice for display, in the mode's physical domain, then min-max to 0..1."""
    arr = slice2d.astype(np.float32)
    m = mode.upper()
    vmax = float(np.nanmax(arr)) if arr.size else 0.0

    if m in {"HA", "IA", "EL"}:
        # 0..255 to −90..90, or assume already degrees
        phys = (arr / 255.0) * 180.0 - 90.0 if vmax > 2.0 else arr
        lo, hi = -90.0, 90.0
    elif m == "AZ":
        # 0..255 to 0..360, or assume already degrees
        phys = (arr / 255.0) * 360.0 if vmax > 2.0 else arr
        lo, hi = 0.0, 360.0
    elif m == "FA":
        # 0..255 to 0..1
        phys = (arr / 255.0) if vmax > 2.0 else arr
        lo, hi = 0.0, 1.0
    else:
        phys = arr
        lo, hi = float(np.nanmin(arr)), float(np.nanmax(arr))

    # normalize for colormap
    norm = (phys - lo) / max(hi - lo, 1e-6)
    return np.clip(norm, 0.0, 1.0)


def plot_label_and_limits(mode: str):
    m = mode.upper()
    if m == "FA":
        return "Fractional Anisotropy", 0.0, 1.0
    if m in {"HA", "IA", "EL"}:
        name = {"HA": "Helical Angle", "IA": "Intrusion Angle", "EL": "Elevation Angle"}[
            m
        ]
        return f"{name} (°)", -90.0, 90.0
    if m == "AZ":
        return "Azimuth (°)", 0.0, 360.0
    return f"{m}", None, None


class Window(QWidget):
    def __init__(
        self,
        output_dir: str,
        mask_path: str | None = "",
        N_slice: int | None = None,
        N_line: int = 5,
        angle_range: float = 20,
        image_mode: str | None = None,
        cmap: Colormap | str | None = None,
    ) -> None:
        super().__init__()

        self.half_point_size = 5
        self.point_size = 10
        self.line_np = None
        self.color_idx = 0
        self.bg_img: QGraphicsPixmapItem | None = None
        self.is_mouse_down = False
        self.start_point: QGraphicsEllipseItem | None = None
        self.end_point: QGraphicsEllipseItem | None = None
        self.start_pos: tuple[float, float] | None = None
        self.end_points: np.ndarray | None = None

        self.N_slice = N_slice or 0
        self.angle_range = angle_range
        self.N_line = N_line
        self.intensity_profiles: list[np.ndarray | None] = []
        self.x_min_lim = None
        self.x_max_lim = None
        self.y_min_lim = None
        self.y_max_lim = None

        self.view = QGraphicsView()
        self.view.setRenderHint(QPainter.Antialiasing)
        self.view.setResizeAnchor(QGraphicsView.AnchorViewCenter)
        self.view.setTransformationAnchor(QGraphicsView.AnchorViewCenter)

        self.mask_path = mask_path
        self.output_path = Path(output_dir)

        # Discover available modes, select default
        self.available_modes = discover_modes(self.output_path)
        if not self.available_modes:
            sys.exit(
                f"No angle or FA folders found in {self.output_path}. Expected one of HA, IA, EL, AZ, FA."
            )
        self.image_mode = (image_mode or self.available_modes[0]).upper()
        if self.image_mode not in self.available_modes:
            sys.exit(
                f"Requested mode {self.image_mode} not found. Available: {self.available_modes}"
            )

        # DataReader for current mode
        self.data_reader = DataReader(self.output_path / self.image_mode)

        # choose default colormap by mode if none provided
        if cmap is None:
            self.cmap = default_cmap_for(self.image_mode)
        elif isinstance(cmap, str):
            self.cmap = plt.get_cmap(cmap)
        else:
            self.cmap = cmap

        if self.N_slice > self.data_reader.shape[0] - 1:
            raise IndexError(
                f"Selected slice {self.N_slice} exceeds number of slices {self.data_reader.shape[0]}"
            )

        self.current_img = self.load_image(self.N_slice)
        self.load_image_to_gui(self.current_img)

        # UI
        menu_bar = QMenuBar(self)
        file_menu = menu_bar.addMenu("File")
        if file_menu is not None:
            act = QAction("Change Slice", self)
            act.triggered.connect(self.change_slice_dialog)
            file_menu.addAction(act)
            act = QAction("Choose image type", self)
            act.triggered.connect(self.change_img_mode)
            file_menu.addAction(act)

        graph_menu = menu_bar.addMenu("Graph")
        if graph_menu is not None:
            a = QAction("Set x axis limits", self)
            a.triggered.connect(self.set_x_lim_dialog)
            graph_menu.addAction(a)
            a = QAction("Set y axis limits", self)
            a.triggered.connect(self.set_y_lim_dialog)
            graph_menu.addAction(a)

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

        # Keep a comfortable default window even for very small input images
        self.setMinimumSize(900, 700)

        self.quit_shortcut = QShortcut(QKeySequence("Ctrl+Q"), self)
        self.quit_shortcut.activated.connect(lambda: quit())

        plot_profile_button.clicked.connect(self.plot_profile)
        save_profile_button.clicked.connect(self.save_profile)

    # ---------- logic ----------

    def update_text(self) -> None:
        angle_range_text = self.input_angle_range.text().strip()
        self.angle_range = (
            float(angle_range_text)
            if angle_range_text and float(angle_range_text) > 0
            else 0.1
        )
        if not angle_range_text or float(angle_range_text) <= 0:
            self.input_angle_range.setText("0.1")

        n_line_text = self.input_N_line.text().strip()
        self.N_line = int(n_line_text) if n_line_text and int(n_line_text) >= 1 else 1
        if not n_line_text or int(n_line_text) < 1:
            self.input_N_line.setText("1")

        self.update_plot()

    def update_plot(self) -> None:
        if not self.start_point or not self.end_point or self.start_pos is None:
            return

        # clear lines
        if hasattr(self, "lines") and self.lines:
            for l in self.lines:
                if l.scene() == self.scene:
                    self.scene.removeItem(l)
        self.lines: list[QGraphicsLineItem] = []

        sx, sy = self.start_pos
        ex, ey = self.end_point.rect().center().x(), self.end_point.rect().center().y()

        self.end_points = find_end_points(
            (sy, sx), (ey, ex), float(self.angle_range), int(self.N_line)
        )

        for i in range(self.N_line):
            line = self.scene.addLine(
                sx,
                sy,
                self.end_points[i][1],
                self.end_points[i][0],
                pen=QPen(QColor("black"), 2),
            )
            if line is not None:
                self.lines.append(line)

        # redraw endpoints
        for item in (self.end_point, self.start_point):
            if item is not None and item.scene() == self.scene:
                self.scene.removeItem(item)

        self.start_point = self.scene.addEllipse(
            sx - self.half_point_size,
            sy - self.half_point_size,
            self.point_size,
            self.point_size,
            pen=QPen(QColor("black")),
            brush=QBrush(QColor("black")),
        )
        self.end_point = self.scene.addEllipse(
            ex - self.half_point_size,
            ey - self.half_point_size,
            self.point_size,
            self.point_size,
            pen=QPen(QColor("black")),
            brush=QBrush(QColor("black")),
        )

    def plot_profile(self) -> None:
        if self.line_np is None or not self.line_np.any():
            print("No line drawn")
            return

        start_point = self.line_np[0:2] * self.bin_factor
        end_point = self.line_np[2:] * self.bin_factor

        label_y, minimum, maximum = plot_label_and_limits(self.image_mode)

        self.intensity_profiles = calculate_intensities(
            self.current_img,
            start_point,
            end_point,
            self.angle_range,
            self.N_line,
            max_value=maximum,
            min_value=minimum,
        )

        plot_intensity(
            self.intensity_profiles,
            label_y=label_y,
            x_max_lim=self.x_max_lim,
            x_min_lim=self.x_min_lim,
            y_max_lim=self.y_max_lim,
            y_min_lim=self.y_min_lim,
        )

    def save_profile(self) -> None:
        if self.line_np is None or not self.line_np.any():
            print("No line drawn")
            return

        start_point = self.line_np[0:2] * self.bin_factor
        end_point = self.line_np[2:] * self.bin_factor

        if not self.intensity_profiles:
            label_y, minimum, maximum = plot_label_and_limits(self.image_mode)
            self.intensity_profiles = calculate_intensities(
                self.current_img,
                start_point,
                end_point,
                self.angle_range,
                self.N_line,
                max_value=maximum,
                min_value=minimum,
            )

        save_path, _ = QFileDialog.getSaveFileName(
            self, "Save Profile", "", "Csv Files (*.csv);;All Files (*)"
        )
        if save_path:
            if not save_path.lower().endswith(".csv"):
                save_path += ".csv"
            save_intensity(self.intensity_profiles, save_path)

    def load_image(self, N_slice: int, bin_factor: int | None = None) -> np.ndarray:
        img = self.data_reader.load_volume(start_index=N_slice, end_index=N_slice + 1)[
            0
        ].astype(np.float32)

        if bin_factor:
            img = block_reduce(img, block_size=(bin_factor, bin_factor), func=np.mean)

        # apply mask before conversion and normalization
        if self.mask_path:
            print(f"\nReading mask from {self.mask_path} ...")
            mreader = DataReader(self.mask_path)
            mask_binf = self.data_reader.shape[0] / mreader.shape[0]
            m = mreader.load_volume(
                start_index=int(self.N_slice / mask_binf),
                end_index=int(self.N_slice / mask_binf) + 1,
            )[0].astype(np.float32)
            m = cv2.resize(
                m, (img.shape[1], img.shape[0]), interpolation=cv2.INTER_LINEAR
            )
            if m.shape != img.shape:
                raise ValueError(
                    f"Mask shape {m.shape} does not match slice shape {img.shape}"
                )
            img[m == 0] = np.nan

        return img

    def load_image_to_gui(self, current_img: np.ndarray) -> None:
        h, w = current_img.shape[:2]

        # choose bin so that display fits below 1000 px each side
        self.bin_factor = 1
        while h // self.bin_factor >= 1000 or w // self.bin_factor >= 1000:
            self.bin_factor *= 2
        print(f"Bin factor: {self.bin_factor}")

        img_bin = block_reduce(
            current_img, block_size=(self.bin_factor, self.bin_factor), func=np.mean
        )

        # convert to physical then normalize for display
        norm01 = convert_slice_for_display(img_bin, self.image_mode)

        # colormap
        cmap = self.cmap or default_cmap_for(self.image_mode)
        rgb = cmap(norm01)
        rgb = (rgb[:, :, :3] * 255).astype(np.uint8)
        gray_color = [128, 128, 128]
        rgb[np.isnan(norm01)] = gray_color

        self.current_img_rgb = rgb
        pixmap = np2pixmap(self.current_img_rgb)

        H, W, _ = self.current_img_rgb.shape
        self.scene = QGraphicsScene(0, 0, W, H)
        self.end_point = None
        self.start_point = None
        self.line = None
        self.lines = []
        self.bg_img = self.scene.addPixmap(pixmap)
        if self.bg_img is not None:
            self.bg_img.setPos(0, 0)
        self.scene.setSceneRect(0, 0, W, H)
        self.view.setScene(self.scene)
        self._fit_view_to_scene()

        self.scene.mousePressEvent = self.mouse_press  # type: ignore
        self.scene.mouseMoveEvent = self.mouse_move  # type: ignore
        self.scene.mouseReleaseEvent = self.mouse_release  # type: ignore

    def _fit_view_to_scene(self) -> None:
        if hasattr(self, "scene") and self.scene is not None:
            self.view.fitInView(self.scene.sceneRect(), Qt.KeepAspectRatio)

    def resizeEvent(self, event) -> None:  # type: ignore[override]
        super().resizeEvent(event)
        self._fit_view_to_scene()

    # ---------- mouse handlers ----------

    def mouse_press(self, event: QGraphicsSceneMouseEvent | None) -> None:
        if event is None:
            return
        x, y = event.scenePos().x(), event.scenePos().y()
        self.is_mouse_down = True
        self.start_pos = (x, y)

        try:
            if self.start_point is not None:
                self.scene.removeItem(self.start_point)
        except Exception:
            pass

        self.start_point = self.scene.addEllipse(
            x - self.half_point_size,
            y - self.half_point_size,
            self.point_size,
            self.point_size,
            pen=QPen(QColor("black")),
            brush=QBrush(QColor("black")),
        )

    def mouse_move(self, event: QGraphicsSceneMouseEvent | None) -> None:
        if not self.is_mouse_down or event is None or self.start_pos is None:
            return

        x, y = event.scenePos().x(), event.scenePos().y()

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

        sx, sy = self.start_pos
        self.end_points = find_end_points(
            (sy, sx), (y, x), float(self.angle_range), int(self.N_line)
        )

        if hasattr(self, "lines") and self.lines:
            for l in self.lines:
                if l.scene() == self.scene:
                    self.scene.removeItem(l)

        self.lines = []
        for i in range(self.N_line):
            line = self.scene.addLine(
                sx,
                sy,
                self.end_points[i][1],
                self.end_points[i][0],
                pen=QPen(QColor("black"), 2),
            )
            if line is not None:
                self.lines.append(line)

    def mouse_release(self, event: QGraphicsSceneMouseEvent | None) -> None:
        if event is None or self.start_pos is None:
            return
        x, y = event.scenePos().x(), event.scenePos().y()
        sx, sy = self.start_pos
        self.is_mouse_down = False
        self.line_np = np.array([sy, sx, y, x], dtype=float)

    # ---------- dialogs ----------

    def change_slice_dialog(self) -> None:
        self.dialog = QWidget()
        self.dialog.setWindowTitle("Change Slice")
        layout = QVBoxLayout()
        self.spin_box = QSpinBox()
        self.spin_box.setMinimum(0)
        self.spin_box.setMaximum(self.data_reader.shape[0] - 1)
        self.spin_box.setValue(self.N_slice)
        layout.addWidget(self.spin_box)
        change_button = QPushButton("Change")
        change_button.clicked.connect(self.change_slice)
        layout.addWidget(change_button)
        self.dialog.setLayout(layout)
        self.dialog.show()

    def change_slice(self) -> None:
        self.N_slice = self.spin_box.value()
        self.current_img = self.load_image(self.N_slice).astype(float)
        self.load_image_to_gui(self.current_img)
        self.dialog.close()

    def set_x_lim_dialog(self) -> None:
        self.dialog = QWidget()
        self.dialog.setWindowTitle("Change X-Axis Limits")
        layout = QVBoxLayout()
        self.x_min_spin_box = QSpinBox()
        self.x_max_spin_box = QSpinBox()
        self.x_min_spin_box.setRange(-10000, 10000)
        self.x_max_spin_box.setRange(-10000, 10000)
        layout.addWidget(QLabel("X Min:"))
        layout.addWidget(self.x_min_spin_box)
        layout.addWidget(QLabel("X Max:"))
        layout.addWidget(self.x_max_spin_box)
        change_button = QPushButton("Set X Limits")
        change_button.clicked.connect(self.set_x_lim)
        layout.addWidget(change_button)
        self.dialog.setLayout(layout)
        self.dialog.show()

    def set_x_lim(self) -> None:
        self.x_min_lim = self.x_min_spin_box.value()
        self.x_max_lim = self.x_max_spin_box.value()
        self.dialog.close()

    def set_y_lim_dialog(self) -> None:
        self.dialog = QWidget()
        self.dialog.setWindowTitle("Change Y-Axis Limits")
        layout = QVBoxLayout()
        self.y_min_spin_box = QSpinBox()
        self.y_max_spin_box = QSpinBox()
        self.y_min_spin_box.setRange(-10000, 10000)
        self.y_max_spin_box.setRange(-10000, 10000)
        layout.addWidget(QLabel("Y Min:"))
        layout.addWidget(self.y_min_spin_box)
        layout.addWidget(QLabel("Y Max:"))
        layout.addWidget(self.y_max_spin_box)
        change_button = QPushButton("Set Y Limits")
        change_button.clicked.connect(self.set_y_lim)
        layout.addWidget(change_button)
        self.dialog.setLayout(layout)
        self.dialog.show()

    def set_y_lim(self) -> None:
        self.y_min_lim = self.y_min_spin_box.value()
        self.y_max_lim = self.y_max_spin_box.value()
        self.dialog.close()

    def change_img_mode(self) -> None:
        self.dialog = QWidget()
        self.dialog.setWindowTitle("Choose Image Mode")
        layout = QVBoxLayout()

        # dynamic radios based on discovered modes
        self._mode_buttons = {}
        group = QButtonGroup(self.dialog)
        for m in self.available_modes:
            rb = QRadioButton(m, self.dialog)
            if m == self.image_mode:
                rb.setChecked(True)
            group.addButton(rb)
            layout.addWidget(rb)
            self._mode_buttons[m] = rb

        confirm_button = QPushButton("Confirm", self.dialog)
        confirm_button.clicked.connect(self._apply_img_mode_selection)
        layout.addWidget(confirm_button)

        self.dialog.setLayout(layout)
        self.dialog.show()

    def _apply_img_mode_selection(self) -> None:
        chosen = next(
            (m for m, rb in self._mode_buttons.items() if rb.isChecked()), None
        )
        if chosen is None:
            return
        self.image_mode = chosen
        print(f"Image mode changed to {self.image_mode}")

        self.data_reader = DataReader(self.output_path / self.image_mode)
        if self.N_slice > self.data_reader.shape[0] - 1:
            self.N_slice = 0

        # update colormap for the new mode
        self.cmap = default_cmap_for(self.image_mode)

        self.current_img = self.load_image(self.N_slice)
        self.load_image_to_gui(self.current_img)
        self.dialog.close()
