"""Logger configuration for the bifidotyper package.

This module sets up logging for the entire package with both file and console output,
configurable log levels, and proper log rotation.
"""

import logging
import sys
import os
from logging.handlers import RotatingFileHandler
from pathlib import Path

class BifidoLogger:
    """Logger class for bifidotyper with configurable output and log rotation."""
    
    def __init__(self, log_file='bifidotyper.log', level=logging.INFO):
        """Initialize the logger with the given settings.
        
        Args:
            log_file (str): Path to the log file. Defaults to 'bifidotyper.log'.
            level (int): Logging level. Defaults to logging.INFO.
        """
        self.log_file = log_file
        self.level = level
        self._setup_logger()
        
    def _setup_logger(self):
        """Configure the logger with both file and console handlers."""
        # Create logger
        self.logger = logging.getLogger('bifidotyper')
        self.logger.setLevel(self.level)
        
        # Remove any existing handlers
        self.logger.handlers.clear()
        
        # Create formatters
        file_formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        console_formatter = logging.Formatter(
            '%(levelname)s - %(message)s'
        )
        
        # File handler with rotation (max 5MB per file, keep 3 backup files)
        file_handler = RotatingFileHandler(
            self.log_file,
            maxBytes=5*1024*1024,  # 5MB
            backupCount=3,
            encoding='utf-8'
        )
        file_handler.setLevel(self.level)
        file_handler.setFormatter(file_formatter)
        
        # Console handler (only INFO and above)
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(console_formatter)
        
        # Add handlers to logger
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)
        
        # Set matplotlib to WARNING level to reduce noise
        logging.getLogger('matplotlib').setLevel(logging.WARNING)
        
    def set_level(self, level):
        """Change the logging level.
        
        Args:
            level (int): New logging level (e.g., logging.DEBUG, logging.INFO)
        """
        self.level = level
        self.logger.setLevel(level)
        for handler in self.logger.handlers:
            handler.setLevel(level)

# Create default logger instance
logger = BifidoLogger().logger

if __name__ == "__main__":
    logger.info("Logger is configured and ready to use.")