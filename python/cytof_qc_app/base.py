"""Shared abstract interface for QC objects."""
from abc import ABC, abstractmethod


class QCObject(ABC):
    def __init__(self, args):
        self.args = args

    @abstractmethod
    def run(self, *data):
        """Run this workflow stage."""

