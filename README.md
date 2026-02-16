# Week Stability Boundary Transfer Tool
Author: Marvin Ahlborn</br>
Date: 2026-02-16

This is a Bi-Circular Restricted Four Body Problem (BCRFBP) simulation tool with two scripts to find Week Stability Boundary (WSB) transfers with ballistic captures to Moon orbits. An example of a transfer and a capture are shown in the following:

<img src="https://github.com/MarvinAhl/wsb-tool/blob/main/plots/transfer.png" alt="A transfer found with the WSB Tool" width="400"/>
<img src="https://github.com/MarvinAhl/wsb-tool/blob/main/plots/capture.png" alt="A ballistic capture found with the WSB Tool" width="300"/>

The BCRFBP equations are taken from Simó et al. [1]. The WSB captures are found using a method similar to Belbruno [2] taken from Luo and Topputo [3]. Additional scripts are available for matching BCRFBP sun angles with real dates and performing a preliminary ground coverage analysis using JPL's SPICE Toolkit [4]. Add SPICE kernels described in the metakernel inside the kernel folder to execute spice related scripts.

[1] C. Simó, G. Gómez, À. Jorba, and J. Masdemont, "The bicircular model near the triangular libration points of the RTBP," From Newton To Chaos, vol. 336, pp. 343-370, 1995, doi: 10.1007/978-1-4899-1085-1_34.</br>
[2] E. Belbruno, "Lunar capture orbits, a method of constructing earth moon trajectories and the lunar GAS mission," in 19th International Electric Propulsion Conference, Colorado Springs, CO, USA: American Institute of Aeronautics and Astronautics, May 1987, doi: 10.2514/6.1987-1054.</br>
[3] Z.-F. Luo and F. Topputo, "Mars orbit insertion via ballistic capture and aerobraking," Astrodyn, vol. 5, no. 2, pp. 167-181, 2021, doi: 10.1007/s42064-020-0095-4.</br>
[4] C. Acton, N. Bachman, B. Semenov, and E. Wright, "A look towards the future in the handling of space science mission geometry," Planetary and Space Science, vol. 150, pp. 9-12, 2018, doi: 10.1016/j.pss.2017.02.013.
