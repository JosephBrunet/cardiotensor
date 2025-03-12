import nox

# Reuse virtual environments to speed up development
# nox.options.reuse_existing_virtualenvs = True

# Define default Python versions for testing
PYTHON_VERSIONS = ["3.10", "3.11", "3.12"]

# Locations to check for linting, formatting, and type checking
PACKAGE_LOCATIONS = ["src", "tests", "noxfile.py", "pyproject.toml"]


@nox.session(python=PYTHON_VERSIONS)
def tests(session):
    """
    Run the test suite using pytest with coverage.
    """
    # Install development dependencies
    session.install(".[dev]")

    # Run pytest without coverage
    session.run("pytest", "-s", *session.posargs)

    # # Run pytest with coverage tracking
    # session.run(
    #     "coverage", "run", "--source=src", "-m", "pytest", "-s", *session.posargs
    # )

    # # Generate coverage reports
    # session.run("coverage", "report", "--show-missing")
    # session.run("coverage", "xml")  # For CI integrations
    # session.run("coverage", "html")  # For local inspection


# @nox.session(python=PYTHON_VERSIONS)
# def lint(session):
#     """
#     Run linters and formatters.
#     """
#     session.install("ruff", "black")
#     session.run("ruff", *PACKAGE_LOCATIONS)
#     session.run("black", "--check", *PACKAGE_LOCATIONS)

# @nox.session(python=PYTHON_VERSIONS)
# def type_check(session):
#     """
#     Run type checking using mypy.
#     """
#     session.install("mypy")
#     session.run("mypy", *PACKAGE_LOCATIONS)
