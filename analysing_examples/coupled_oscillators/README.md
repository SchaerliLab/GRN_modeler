# Intracellularly Coupled Oscillators for Synthetic Biology

This repository contains the simulation scripts and models accompanying our manuscript:

**"Intracellularly Coupled Oscillators for Synthetic Biology"**  
submitted to *Nature Communications*.

The code includes deterministic, stochastic, and hybrid models of coupled genetic oscillators, exploring synchronization, frequency and amplitude modulation, resonance, and bifurcation behavior.

---

## Requirements

- **MATLAB** with the **SimBiology Toolbox** (all models rely on it).
- **COPASI**, **Python**, and **BasiCO** (for Lyapunov exponent calculations).
- **C++ compiler** (for accelerating models and stochastic simulations).

---

## Repository Structure

### Deterministic Simulations (MATLAB scripts)

- `computer_CRISPRlator.m`  
  Deterministic biocomputers, phase locking through dCas system and nonlinear response function.

- `coupled_synchronisation.m`  
  Weakly coupled repressilators, coupled through a common protease system.

- `deeply_coupled_prot.m`  
  Two repressilators coupled through a common edge with varying protease concentrations; bifurcation diagrams and Poincaré maps.

- `deeply_coupled_tirv.m`  
  Two repressilators coupled through a common edge, one node strength varied with an extra protease system.

- `deeply_coupled_tirv_Elowitz.m`  
  Two repressilators coupled through a common edge with the Elowitz model and varying promoter strengths.

- `deeply_coupled_tirv_Elowitz_inducer.m`  
  Two repressilators coupled through a common edge, repression strength modulated with an external inducer.

- `double_repressilator_T.m`  
  Calculates oscillation periods of deeply coupled repressilators with different node numbers.

- `frequency_modulation.m`  
  Repressilator coupled with a Goodwin oscillator, demonstrating frequency modulation.

- `goodwin.m`  
  Goodwin oscillator oscillatory region and bifurcation diagram.

- `independent2.m`  
  CRISPRlator and Repressilator repressing a node; frequency modulation.

- `independent3.m`  
  CRISPRlator and Stricker oscillator repressing a node; frequency modulation.

- `independent_oscillators.m`  
  Two repressilators repressing a common node; beat phenomenon.

- `modulation_independent.m`  
  Circuit for independent amplitude modulation (Tomazou design).

- `repressilator_goodwin.m`  
  Repressilator with a Goodwin oscillator input; Lyapunov exponent, bifurcation diagram, Poincaré map.

- `repressilator_prot.m`  
  Repressilator with changing protease concentration; Lyapunov exponent, bifurcation diagram, Poincaré map.

- `repressilator_repressilator.m`  
  Unidirectionally coupled repressilators; Lyapunov exponent, bifurcation diagram, Poincaré map.

- `repressilator_T.m`  
  Repressilators with different node numbers; oscillation period calculation.

- `resonance.m`  
  Unidirectionally coupled repressilators showing resonance.

- `togglelator.m`  
  Oscillatory behavior of the toggle switch; bifurcation diagram compared with the Goodwin oscillator.

- `togglelator_phasespace.m`  
  Phase space for the toggle switch.

- `Tomazou_repressilator_T.m`  
  Oscillation period of repressilators as a function of node number (Tomazou model).

---

### Stochastic Simulations (`computer_stoch/`)

These simulations require a C++ backend for efficiency.

- `computer_CRISPRlator_sctoch_burden.m`  
  Effect of metabolic burden in coupled CRISPRlators with a common dCas system.

- `computer_CRISPRlator_sctoch_repress.m`  
  Effect of explicit interactions in coupled CRISPRlators under varying noise strength.

---

## Citation

If you use this code in your research, please cite our paper:  
**"Intracellularly Coupled Oscillators for Synthetic Biology"**, *Nature Communications*, (202X).

---

## Contact

For questions or issues, please open a GitHub issue or contact:  
*Gabor Hollo* (hollo88@gmail.com)
