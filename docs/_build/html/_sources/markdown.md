# Setup

Set up your environment on Linux. Alternatively, use WSL2 on Windows or macOS.
```{note}
The scripts were developed and tested on Windows using WSL2.
```
**info:** Here is detailed the setup on Windows using WSL2.

## What is WSL2?
Windows Subsystem for Linux allows developers to run a Linux environment on Windows without the need for a separate **virtual machine** or **dual booting**. 
WSL gives you access to a full Ubuntu terminal environment and allows you to work with different Linux distributions.
```{note}
To set up WSL2: [Install Linux on Windows with WSL](https://learn.microsoft.com/en-us/windows/wsl/install)
```

## Terminals and Distributions

You can either use PowerShell(pre-installed on Windows) and Ubuntu terminal to
access the Linux environment. Installing [Ubuntu 22.04.2 LTS](https://apps.microsoft.com/store/detail/ubuntu-22042-lts/9PN20MSR04DW?hl=en-gb&gl=gb&rtc=1) will provide you with both the terminal and Linux distribution.

```{note}
Three terminals were used while developing the scripts: PowerShell, Ubuntu and Anaconda.
Powershell: used to install dependencies that failed to install on Anaconda. </br>
Ubuntu: access the development environment for the scripts. </br>
Anaconda: primary point for installing the dependencies for developing the scripts.
```

## Dependencies
| Name                | Version  | Build              | Channel      |
| ------------------- | -------- | ------------------ | ------------ |
| <span style="color:blue"/>**ase**             | 3.22.1   | pyhd8ed1ab_1       | conda-forge  |
| <span style="color:blue"/>**click**           | 8.1.3    | win_pyhd8ed1ab_2   | conda-forge  |
| <span style="color:blue"/>**matplotlib-base** | 3.6.2    | py310h51140c5_0    | conda-forge  |
| <span style="color:blue"/>**numpy**           | 1.23.5   | py310h4a8f9c9_0    | conda-forge  |
| <span style="color:blue"/>**pandas**          | 1.5.2    | pypi_0             | pypi         | 
| <span style="color:blue"/>**pip**             | 22.3.1   | pyhd8ed1ab_0       | conda-forge  |
| <span style="color:blue"/>**python**          | 3.10.6   | hcf16a7b_0_cpython | conda-forge  |
| <span style="color:blue"/>**seaborn**         | 0.12.2   | pypi_0             | pypi         |
| <span style="color:blue"/>**shiftml2**        | 1.0.0    | pypi_0             | pypi         |
| <span style="color:blue"/>**soprano**         | 0.8.13   | pypi_0             | pypi         |
| <span style="color:blue"/>**wheel**           | 0.38.4   | pyhd8ed1ab_0       | conda-forge  |

## Learn more

This is just a simple starter to get you started.
You can learn a lot more at [jupyterbook.org](https://jupyterbook.org).
