# Plot_RNA2



# Installation Instructions

To install the software, follow these steps:

1. **Create a Conda Environment:**
   Open your terminal or command prompt and run the following command to create a new Conda environment named as you want, here `Plot_RNA` with Python version 3.13.2:
   ```bash
   conda create --name Plot_RNA python=3.13.2
   ```

2. **Activate the Conda Environment:**
   Activate the newly created environment by running:
   ```bash
   conda activate Plot_RNA
   ```

3. **Create a Directory for the program:**
   Create a directory named `RNA_plots` to save the program:
   ```bash
   mkdir RNA_plots
   ```

4. **Run the Installation Script:**
   Execute the following command to download and run the installation script:
   ```bash
   python <(curl -s https://raw.githubusercontent.com/louismeuret/Plot_RNA2/refs/heads/main/install_script.py) --repo https://github.com/louismeuret/Plot_RNA2 --install-path .
   ```

# Launching the Server

To launch the server, use one of the following commands based on your operating system:

- **For Linux/MacOS:**
  ```bash
  ./run_program.sh
  ```

- **For Windows:**
  ```bash
  run_program.bat
  ```
