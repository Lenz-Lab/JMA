# Dynamic Joint Measurement Analysis
![Distance_TTA_vs_Control](https://user-images.githubusercontent.com/69816397/211891169-d905490a-9692-4162-b5f7-3c27e6a9c24a.gif)
![TTA_Representative](https://github.com/Lenz-Lab/JMA/assets/69816397/84325e31-6759-422e-a0de-fb249d2b93e1.gif)
## Description
The Joint Measurement Analysis (JMA) toolbox is a set of MATLAB scripts. It measures joint space distance and congruence between two bones at statistical shape model correspondence particles, either statically or through a dynamic activity. It then compares groups and shows the results on the mean bone surface as images and videos.

For detailed, step-by-step instructions, see the standard operating procedure: [SOP_Joint_Measurement_Analysis.docx](SOP_Joint_Measurement_Analysis.docx).

## Workflow
![JMA_Workflow](https://github.com/user-attachments/assets/61e2161b-e3a9-4342-91f2-40ce85ed5e20)

Run the scripts in order from the repository folder. Each one prompts for its inputs through dialogs, and adds [Scripts/](Scripts/) to the MATLAB path itself.

| Script | What it does | Main output |
| --- | --- | --- |
| [Scripts/JMA_00_DSX_PreProcess.m](Scripts/JMA_00_DSX_PreProcess.m) | Optional. Converts DSX transform exports into the per-bone kinematics `.txt` files. | Kinematics `.txt` files |
| [JMA_01_Kinematics_to_SSM.m](JMA_01_Kinematics_to_SSM.m) | Applies the kinematics to each subject's bones and calculates distance and congruence at every correspondence particle, frame by frame. | `<Subject>\Data_<Bone1>_<Bone2>_<Subject>.mat`, `Outputs\JMA_01_Outputs\Data_<Bone1>_<Bone2>.mat` |
| [JMA_01a_Import_User_Data.m](JMA_01a_Import_User_Data.m) | Optional. Adds your own per-particle data (e.g. FEA, cortical thickness) from `.xlsx`/`.csv` to the JMA_01 outputs. | Updated JMA_01 `.mat` files |
| [JMA_02_Data_Process_and_Normalize.m](JMA_02_Data_Process_and_Normalize.m) | Normalizes every subject to percent of stance and pools the data by group. | `Outputs\JMA_02_Outputs\Normalized_Data_<Bone1>_<Bone2>*.mat` |
| [JMA_03_Statistical_Analyses.m](JMA_03_Statistical_Analyses.m) | Group statistics and visualization on the mean bone (see below). | `Results\` folder with `.tif`, `.mp4` and `.xlsx` files |
| [JMA_04_Dynamic_Visualization.m](JMA_04_Dynamic_Visualization.m) | Shows group or individual results on the bones while they move. | `.mp4` |

## Requirements
- MATLAB R2019b or newer (developed on R2023a)
- Statistics and Machine Learning Toolbox
- Simulink 3D Animation (`vrrotvec`, used to orient the glyphs in the result figures)
- Parallel Computing Toolbox (JMA_01, and regional percentages in JMA_03's SPM mode)
- [spm1d](https://spm1d.org/install/InstallationMatlab.html) on the MATLAB path (only for Statistical Parametric Mapping in JMA_03)
- Correspondence particles and mean shapes from [ShapeWorks](https://sfi.utah.edu/software/shapeworks/)

## Data Folder Structure
Spelling of bone and group names must match across every file name.

```
<Main Directory>
├── <Group_A>
│   ├── <Subject_01>
│   │   ├── <Name>.local.particles   (from ShapeWorks)
│   │   ├── <Bone_Name_01>.stl       (bone the data is mapped to)
│   │   ├── <Bone_Name_02>.stl       (opposing bone)
│   │   ├── <Name>.xlsx              (gait events: first tracked, heel-strike, toe-off, last tracked frame)
│   │   ├── <Bone_Name_01>.txt       (4x4 transform per frame, one line per frame)
│   │   └── <Bone_Name_02>.txt
│   └── <Subject_...>
├── <Group_...>
├── Mean_Models                      (needed for JMA_03)
│   ├── <Group>_<Bone_Name_01>.stl
│   └── <Group>_<Bone_Name_01>.particles
└── Outputs                          (created by the scripts)
```

JMA_03 finds each group's mean model by splitting the file name on `_` and matching both the group name and the bone name. For example, `Shod_Calcaneus.stl` works for group `Shod`, but `Insole_32_Calcaneus.stl` does **not** match group `Insole32`.

## JMA_03 Analysis Modes
| Mode | Groups | Description |
| --- | --- | --- |
| 1. Statistical Analyses | 2 or more | Per-particle, per-frame tests. **Two-Sample:** t-test and Wilcoxon rank sum, or with *Paired Data* checked, a paired t-test and signed rank test. Runs every pair of groups in both directions. **Multi-Group:** one-way ANOVA and Kruskal-Wallis with post-hoc comparison of the pair you select. |
| 2. Statistical Parametric Mapping | exactly 2 | SPM (spm1d) across stance at each particle. Can optionally report the percentage of significant particles within `.stl` regions of the mean bone. Dynamic data only. |
| 3. Visualization Only: Group | 1 | Group mean at each particle, no statistics. |
| 4. Visualization Only: Individual | 1 group, chosen subjects | Each subject's own values on their own bone. Can optionally align the bone to the particles with ICP. |
| 5. Error | 1 | Absolute percent error of one measure compared with a "ground truth" measure. |

Settings chosen in the first dialog:
- **Frame Rate**: frame rate of the output `.mp4` (dynamic data only).
- **Combine statistical analyses** (mode 1): if checked, a particle is significant when either the parametric or the nonparametric test is. If unchecked, only the test that matches the normality result is used.
- **Minimum percentage of participants** (modes 1, 3, 5): a particle is only tested or drawn if at least this percentage of each group's subjects has data there.

The limits dialog then sets the colorbar range for each measure. It also sets a distance cutoff: particles whose mean distance falls outside the cutoff are left out.

The first figure opens a figure-settings editor for view, colors and glyphs. These settings can be saved to `Outputs\JMA_03_Outputs` and loaded on later runs.

**Outputs**, written to `<Main Directory>\Results\`:
- `<Test>_<Measure>_<Bones>\...\*.tif`: one image per frame (percent of stance)
- `<Test>_<Measure>_<Bones>\*.mp4`: frames stitched into a video (dynamic data)
- `<Measure>_Distributions_<Bones>_FullJoint.xlsx` and `<Measure>_EffectSize_<Bones>_FullJoint.xlsx` (modes 1–2)
- `SPM_Percentages\` (mode 2 with regions)

## Citation
If you use this toolbox, please cite it as described in [CITATION.cff](CITATION.cff).

## License
See [LICENSE](LICENSE). Licenses for the third-party functions in [Scripts/](Scripts/) are in [Scripts/Script_Licenses/](Scripts/Script_Licenses/).
