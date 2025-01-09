import nox

# Define default Python versions for testing
PYTHON_VERSIONS = ["3.10", "3.11", "3.12"]

# Locations to check for linting, formatting, and type checking
PACKAGE_LOCATIONS = ["src", "tests", "noxfile.py", "pyproject.toml"]


@nox.session(python=PYTHON_VERSIONS)
def tests(session):
    """
    Run the test suite using pytest.
    """
    # Install dependencies
    session.install("pytest", "pytest-cov")
    # Run tests with coverage
    session.run("pytest", "--cov=src", "--cov-report=term-missing", *session.posargs)


@nox.session
def lint(session):
    """
    Lint the codebase using flake8.
    """
    session.install("flake8")
    session.run("flake8", *PACKAGE_LOCATIONS)


@nox.session
def format_code(session):
    """
    Format the codebase using black.
    """
    session.install("black")
    session.run("black", "--check", *PACKAGE_LOCATIONS)


@nox.session
def type_check(session):
    """
    Run type checking using mypy.
    """
    session.install("mypy")
    session.run("mypy", *PACKAGE_LOCATIONS)


@nox.session
def safety_check(session):
    """
    Check for vulnerabilities in dependencies using safety.
    """
    session.install("safety")
    session.run("safety", "check", "--full-report")


@nox.session
def docs(session):
    """
    Build the documentation using Sphinx.
    """
    session.install("sphinx", "sphinx-rtd-theme")
    session.run("sphinx-build", "docs", "docs/_build")


@nox.session
def format_fix(session):
    """
    Auto-format the codebase using black.
    """
    session.install("black")
    session.run("black", *PACKAGE_LOCATIONS)
