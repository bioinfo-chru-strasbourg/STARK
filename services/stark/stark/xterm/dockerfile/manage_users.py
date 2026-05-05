from passlib.context import CryptContext
import json

# --- Configuration ---
USERS_FILE = "users.json"
pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")

def load_users():
    """Loads users from the JSON file."""
    try:
        with open(USERS_FILE, "r") as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return {}

def save_users(users_data):
    """Saves users to the JSON file."""
    with open(USERS_FILE, "w") as f:
        json.dump(users_data, f, indent=4)

def hash_password(password: str):
    """Hashes a password."""
    return pwd_context.hash(password)

def add_user(username, password):
    """Adds a new user with a hashed password."""
    users = load_users()
    if username in users:
        print(f"Error: User '{username}' already exists.")
        return
    users[username] = hash_password(password)
    save_users(users)
    print(f"User '{username}' added successfully.")

def remove_user(username):
    """Removes a user."""
    users = load_users()
    if username not in users:
        print(f"Error: User '{username}' not found.")
        return
    del users[username]
    save_users(users)
    print(f"User '{username}' removed successfully.")

def list_users():
    """Lists all users."""
    users = load_users()
    if not users:
        print("No users found.")
        return
    print("Users:")
    for username in users:
        print(f"- {username}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Manage users for the STARK Web Terminal.")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Add user command
    parser_add = subparsers.add_parser("add", help="Add a new user")
    parser_add.add_argument("username", help="The username to add")
    parser_add.add_argument("password", help="The password for the new user")

    # Remove user command
    parser_remove = subparsers.add_parser("remove", help="Remove a user")
    parser_remove.add_argument("username", help="The username to remove")

    # List users command
    subparsers.add_parser("list", help="List all users")

    args = parser.parse_args()

    if args.command == "add":
        add_user(args.username, args.password)
    elif args.command == "remove":
        remove_user(args.username)
    elif args.command == "list":
        list_users()
    else:
        parser.print_help()
