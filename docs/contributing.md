# Contributing

Thank you for considering contributing to the **cardiotensor** project! We welcome contributions of all kinds, including bug fixes, feature suggestions, documentation improvements, and more. Follow the guidelines below to ensure a smooth contribution process.

## How to Contribute

1. **Fork the Repository**
    - Go to the [GitHub repository](https://github.com/JosephBrunet/cardiotensor) and click the "Fork" button.
    - Clone your fork locally:
        ```bash
        git clone https://github.com/<your-username>/cardiotensor.git
        cd cardiotensor
        ```

2. **Create a Branch**
    - Create a new branch for your feature or bug fix:
        ```bash
        git checkout -b my-feature-branch
        ```

3. **Install Dependencies**
    - Ensure you have all required dependencies installed. You can use the provided development dependencies:
        ```bash
        pip install -e .[dev]
        ```

4. **Make Your Changes**
    - Make your changes to the codebase, documentation, or both.
    - Follow [PEP8](https://peps.python.org/pep-0008/) standards for Python code.

5. **Run Tests**
    - Ensure all tests pass and that your contribution does not introduce any issues:
        ```bash
        pytest
        ```
    - If you add new features, include corresponding tests.

6. **Commit Your Changes**
    - Write clear and concise commit messages:
        ```bash
        git add .
        git commit -m "Add feature: XYZ"
        ```

7. **Push and Create a Pull Request**
    - Push your branch to your fork:
        ```bash
        git push origin my-feature-branch
        ```
    - Go to the original repository and create a Pull Request (PR).
    - Provide a clear description of your changes and why they are necessary.

## Coding Standards

To maintain a consistent and readable codebase, please adhere to the following guidelines:

- **PEP 8**: Follow [PEP8](https://peps.python.org/pep-0008/) coding style standards for Python code.
- **Type Annotations**: Include type hints for function arguments and return values.
- **Docstrings**: Provide clear docstrings for functions and classes. We recommend using the Google docstring style.

## Guidelines

- **Be Respectful**: Treat all contributors with respect and professionalism.
- **Follow Standards**: Ensure your contributions align with the project's coding and documentation standards.
- **Keep It Simple**: Focus on making the project easier to use and maintain.
- **Documentation**: Update documentation if your changes affect the usage or functionality of the project.

## Reporting Issues

If you encounter a bug or have a suggestion, please [open an issue](https://github.com/JosephBrunet/cardiotensor/issues). Include as much detail as possible:

- Steps to reproduce the issue.
- Expected and actual results.
- Your environment (e.g., Python version, OS).

## Contact

For questions or further assistance, feel free to reach out to the project maintainers via the repository's [Discussions](https://github.com/JosephBrunet/cardiotensor/discussions) or email listed in the repository.

We appreciate your contributions and support in making cardiotensor better!

