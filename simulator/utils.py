import json
import os

# Get the absolute path of the config file inside the package
CONFIG_PATH = os.path.join(os.path.dirname(__file__), "config.json")

def load_config():
    """Load JSON configuration file."""
    with open(CONFIG_PATH, "r") as f:
        return json.load(f)
