# 2dIR
1d- and 2d-IR calculation in semi-classical line shape theory.
Build in C++ with Python/Cython interface. 

## useage detail
See Readme in source/ and the test.py for example.

## theory and computational tricks
See PDF file in note/

The basic formula for three-order non-linear spectroscopy is written as:
<img src="https://user-images.githubusercontent.com/38717633/235039449-cb91b494-a934-4aa0-88ee-5204f852042b.png" width=80% height=80%>

Where $R_r$ and $R_{nr}$ means rephasing and non-rephasing response function:

<img src="https://user-images.githubusercontent.com/38717633/235039395-11ffb663-9961-43a3-a659-fb6f2d4837f9.png" width=50% height=50%>


In semi-classical line shape theory, as an example, $R_3$ can be written as:

<img src="https://user-images.githubusercontent.com/38717633/235039983-c4b7e9bd-a2ec-4e23-8dcc-bb07013b9118.png" width=50% height=50%>
<img src="https://user-images.githubusercontent.com/38717633/235040624-1b58c511-34c8-46d4-99f7-4f1214a8e05b.png" width=30% height=30%>

where $\mu_{ij}$, $\omega_{ij}$ is the corresponding transition dipole and transition frequencies, respectively, of i,j states of a given vibration mode. <...> is classical ensemble.

The main feature of this program is calculating $I$ with $\omega_{ij}(t)$ and $\mu_{ij}(t)$ as input.
