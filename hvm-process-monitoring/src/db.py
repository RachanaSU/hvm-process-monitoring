"""Database access.

Uses SQLAlchemy when it is installed (the normal case: pip install -r
requirements.txt). Falls back to Python's built-in sqlite3 driver so the
pipeline still runs on a machine without SQLAlchemy. pandas accepts either.
"""
import sqlite3

from .config import DB_PATH

try:
    from sqlalchemy import create_engine
    HAVE_SQLALCHEMY = True
except ImportError:  # pragma: no cover
    HAVE_SQLALCHEMY = False


def get_connection(path=DB_PATH):
    path.parent.mkdir(parents=True, exist_ok=True)
    if HAVE_SQLALCHEMY:
        return create_engine(f"sqlite:///{path}")
    return sqlite3.connect(path)


def driver_name():
    return "SQLAlchemy" if HAVE_SQLALCHEMY else "sqlite3 (SQLAlchemy not installed)"
