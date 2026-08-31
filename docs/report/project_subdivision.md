# AMSC Project Proposal

The project will consist in the extension of LBM algorithm developed as hands-on.
The proposal is the following (considering that a better collision operator than BGK should be implemented):

- extending the LBM simulation to the 3D case;
- parallelizing the LBM algorithm using CUDA;
- simulating the poiseuille law;
- study of turbulence in wing section;

## Subdivision of Tasks

### The team

- Chiara Nonino;
- Alessandro Frisone;
- Corrado Sciancalepore;
- Gabriele Santandrea.

### Tasks and possible division

Structure:
`- Task (assignees) (branch-name)`

-- Considering Implementing on Cluster --

- Poiseuille flow simulation (Frisone, Santandrea + 1) (sim-poiseuille);
- Collision operator:
    - BGK (Already Implemented);
    - TRT (Sciancalepore) (lbm-trt);
    - ? (Nonino) (lbm-?);
- CUDA implementation of LBM 2D (Nonino + 1) (CUDA-2D);
- 3D extension (Santandrea, Sciancalepore) (lbm-3D);
- CUDA implementation of LBM 3D (Frisone + 1) (CUDA-3D);
- Study of turbulence in wing section (Frisone, Santandrea, Sciancalepore, Nonino) (sim-wing-sec);

--- Hybrid Implementation (CUDA + OpenMP) ---
