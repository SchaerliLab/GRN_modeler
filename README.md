# GRN_Modeler

Welcome to **GRN_Modeler**, an application built on MATLAB SimBiology for constructing and simulating gene regulatory networks!

---

## Reference

**GRN_modeler: An Intuitive Tool for Constructing and Evaluating Gene Regulatory Networks and its Applications to Oscillators and a Light Biosensor**  
DOI: [10.1101/2024.12.18.629005 ](https://www.biorxiv.org/content/10.1101/2024.12.18.629005v2)

---

## Installation

To use the application, follow these steps:

1. **Install MATLAB**  
   Download and install [MATLAB](https://uk.mathworks.com/products/matlab.html) (version 2022b or later) along with the [SimBiology Toolbox](https://uk.mathworks.com/products/simbiology.html).

2. **Clone the Repository**  
   Use the following command in your terminal:  
   ```bash
   git clone https://github.com/hollo88/GRN.git
   ```  
   Alternatively, clone the repository directly through MATLAB by following [these instructions](https://uk.mathworks.com/help/simulink/ug/clone-git-repository.html), or download and extract the files manually.

3. **Run the Application**  
   Start MATLAB and type the following command in the command line:  
   ```matlab
   GRN_modeler
   ```

### Optional: Using COPASI for Solvers

If you plan to use a COPASI solver, additional installations are required:

- Install COPASI: [Download](https://copasi.org/Download/) and [Install](https://copasi.org/Support/Installation/), then add it to your system PATH:
  - [Windows Instructions](https://www.architectryan.com/2018/03/17/add-to-the-path-on-windows-10/)
  - [Mac Instructions](https://www.architectryan.com/2012/10/02/add-to-the-path-on-mac-os-x-mountain-lion/)
  - [Linux Instructions](https://linuxize.com/post/how-to-add-directory-to-path-in-linux/)

- Install a supported version of [Python](https://www.python.org/downloads/) and [pip](https://pip.pypa.io/en/stable/installation/).

- Install the BASICO package to connect MATLAB with COPASI:
  ```bash
  pip install copasi-basico
  ```

---

## Manual

A brief overview of the command-line functionalities is available in the Supplementary Information (SI) section of our [manuscript](http://t1.gstatic.com/licensed-image?q=tbn:ANd9GcS8piVHY8Mf7p62OkwcgqmZjlWI8WS58uafqox6EPKAt6OE8IWu3Aow6W-iQHJJLJ09).  
The SI also includes a video tutorial demonstrating how to use the graphical user interface (GUI).

---

## Examples

The repository contains several [example files](/examples) to help you get started. To load an example, navigate to **Input/Output > Load** within the application and select a file from the `examples` folder.  
Simulations presented in the [manuscript](http://t1.gstatic.com/licensed-image?q=tbn:ANd9GcS8piVHY8Mf7p62OkwcgqmZjlWI8WS58uafqox6EPKAt6OE8IWu3Aow6W-iQHJJLJ09) can be found in the `analysing_examples/case_studies1` folder.

---

## Troubleshooting

1. **Accelerating MATLAB Simulations**  
   Ensure a supported C/C++ compiler is installed. A list of supported compilers is available [here](https://uk.mathworks.com/support/requirements/supported-compilers-linux.html).

2. **Python Compatibility**  
   Check MATLAB's supported Python versions [here](https://uk.mathworks.com/support/requirements/python-compatibility.html).  
   To locate your installed Python version and path:
   - **Windows**:  
     ```bash
     where python3
     ```
   - **Mac/Linux**:  
     ```bash
     which python3
     ```
   Check the version with:
   ```bash
   python3 --version
   ```
   In MATLAB, use the following command to specify the Python version:
   ```matlab
   pyenv('Version', 'python_path');
   ```
   Replace `python_path` with the Python installation path (e.g., `'/usr/bin/python3'` on Linux or `'C:\Users\UserName\AppData\Local\Programs\Python\Python311\python'` on Windows).  
   If needed, switch to **OutOfProcess** execution mode:
   ```matlab
   pyenv(ExecutionMode="OutOfProcess");
   ```

3. **BASICO Module Not Found**  
   Ensure BASICO is installed with the correct pip version. Use the `python_path` identified earlier:
   ```bash
   "python_path" -m pip install copasi-basico
   ```

   Replace `python_path` with the Python installation path.

---

With these steps, you're ready to explore the exciting world of gene regulatory network modeling!
