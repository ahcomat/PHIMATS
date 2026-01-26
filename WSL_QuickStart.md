
## PHIMATS: 15-Minute Environment Setup for Windows Users

Follow these steps in order to prepare your Windows machine for `PHIMATS` and other high-performance scientific computing applications.

### 1️⃣ Enable the Linux Engine (WSL2)

Follow the instructions [here](https://learn.microsoft.com/en-us/windows/wsl/install)

1. Open **PowerShell** as Administrator (Right-click Start button).
2. Type this command and hit Enter:
```powershell
wsl --install
```
3. **Restart your computer.**
4. After restart, an Ubuntu window will pop up. Create your **username** and **password**.

**Recommended terminal** [Windows Terminal app](https://apps.microsoft.com/detail/9n0dx20hk701?hl=en-US&gl=FI)

### 2️⃣ On Windows, install the following

Download and install these apps like any other Windows software:

* **Docker Desktop:** [Download here](https://www.docker.com/products/docker-desktop/)
  * Inside Docker Desktop Settings, go to `Resources > WSL Integration` and ensure "Ubuntu" is turned **ON**. 
  * Instructions for integration with WSL2 [here](https://docs.docker.com/desktop/features/wsl/)
* **Recommended IDE: VS Code** [Download here](https://code.visualstudio.com/)

---

### 3️⃣ Paraview for visualization

**ParaView:** [Download here](https://www.paraview.org/download/)

* If you use Windows installation, recommended version is **5.9.1** to avoid incompatibilities with `.xdmf` file paths in `WSL`.
* For newer versions, a **native Linux** setup is recommended. 

---

### 4️⃣ Setup VS Code for Linux

1. Open **VS Code**.
2. Click the **Extensions** icon (the 4 squares on the left).
3. Search for and install: **"WSL"** and **"Dev Containers"**.
4. In the bottom-left corner, click the **"><" icon** and select **"Connect to WSL"**.

**Optional** install the [`C/C++` extension pack](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools-extension-pack) and [H5Web extension](https://marketplace.visualstudio.com/items?itemName=h5web.vscode-h5web) of inspecting `.hdf5` files.

---

### 5️⃣ Get PHIMATS and Docker Ready

1. In a terminal, **Clone the code:**
```bash
git clone https://github.com/ahcomat/PHIMATS.git
cd PHIMATS
```
2. **Pull the pre-built Engine:**
(Make sure Docker Desktop is running in the background!)
```bash
docker pull ahcomat/phimats_dep:latest
```

---

### ✅ Success!

Your machine is now ready. Whenever you want to work on `PHIMATS`, just open your folder in VS Code, open the terminal, and run your `docker run` command from the `INSTALL.md` instructions.