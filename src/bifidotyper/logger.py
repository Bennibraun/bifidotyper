import logging
import sys

# Configure the logger
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('bifidotyper.log'),
        # logging.StreamHandler(sys.stdout),
    ]
)

# Create a logger instance
logger = logging.getLogger('bifidotyper')

logging.getLogger('matplotlib').setLevel(logging.WARNING)

# Example usage
if __name__ == "__main__":
    logger.info("Logger is configured and ready to use.")