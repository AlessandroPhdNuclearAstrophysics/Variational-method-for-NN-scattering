#! /bin/bash
# This script checks if the required dependencies are installed.
REQUIRED_DEPENDENCIES=(
    "git"
    "python3"
    "diff-numerics"
    "grace"
)

ask_install() {
    echo "Do you want to install it? (y/n)"
    read -r answer
    if [[ "$answer" == "y" || "$answer" == "Y" ]]; then
        case "$1" in
            "git")
                sudo apt-get install git
                ;;
            "python3")
                sudo apt-get install python3
                ;;
            "diff-numerics")
                sudo apt-get install diff-numerics
                ;;
            "grace")
                sudo apt-get install grace
                ;;
            *)  
                echo "Unknown dependency: $1"
                ;;
        esac
    else
        echo "Skipping installation of $1."
    fi
}

check_dependencies() {
    for dependency in "${REQUIRED_DEPENDENCIES[@]}"; do
        if ! command -v "$dependency" &> /dev/null; then
            echo "Error: $dependency is not installed."
            ask_install "$dependency"
            return 1
        fi
    done
    echo "All required dependencies are installed."
    return 0
}

check_dependencies