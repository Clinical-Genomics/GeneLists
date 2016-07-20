import re

def cleanup_description(description):
    """Remove the comment in the description and clean up invalid characters: ,:;|>

    Args:
        description (str): text to clean up

    Returns: str or None

    """
    if description:
        description = description.strip()
        description = re.sub(r'\[.*\]', '', description)
        description = re.sub(r'[,:;>| ]', '_', description)
        if description.endswith('_'):
            description = description[:-1]
        return description
    return ''
