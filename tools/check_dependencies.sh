#! /bin/bash
# This script checks if the required dependencies are installed.
REQUIRED_DEPENDENCIES=(
    "git"
    "python3"
    "pip"
    "diff-numerics"
    "grace"
    "fortran-module-interface-generator"
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
            "pip")
                sudo apt-get install python3-pip --fix-missing
                ;;
            "diff-numerics")
                git clone https://github.com/AlessandroPhdNuclearAstrophysics/diff-numerics
                cd diff-numerics
                make -j
                sudo make install
                cd ..
                rm -rf diff-numerics
                ;;
            "grace")
                sudo apt-get install grace
                ;;
            "fortran-module-interface-generator")
                git clone https://github.com/AlessandroPhdNuclearAstrophysics/fortran-module-public-interface-generator
                cd fortran-module-public-interface-generator
                pip install -r requirements.txt
                python3 setup.py install -e
                cd ..
                rm -rf fortran-module-public-interface-generator
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
    ALL_DEPENDENCIES_INSTALLED=true
    for dependency in "${REQUIRED_DEPENDENCIES[@]}"; do
        if ! command -v "$dependency" &> /dev/null; then
            echo "Error: $dependency is not installed."
            ask_install "$dependency"
            ALL_DEPENDENCIES_INSTALLED=false
        fi
    done
    if [ "$ALL_DEPENDENCIES_INSTALLED" = true ]; then
        echo "All required dependencies are installed."
    fi
    return 0
}

check_dependencies