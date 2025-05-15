CS140 Mechanics Final Project


Contributors:
* Ruben Sargsyan - #A09220405
* Vahe Jraghatspanyan - #A09220179
* Kevin Markalian - #A09220549


Format:
* The project will be done by a group, since individual submissions will be suboptimal and will take longer to implement. We decided to combine our mechanical knowledge and programming skills to get this project done.
* Vahe had technical issues regarding the Eclipse IDE and Kevin was heavily loaded with exams, so working in a group seemed more doable.


Task distributions:
* Ruben Sargsyan and Vahe Jraghatspanyan will work on the coding part.
* Kevin Markalian will write the summary report about the tasks and will help out with theoretical parts of the project.


________________


________________


Final Report




Overview
This project explored the formulation and implementation of 2D packing and floor-planning problems using particle-based simulations solved by implementing classical mechanics formulas. By representing geometric items (circles, squares, rectangles) as interacting particles with attractive and repulsive forces, we simulated natural convergence behaviors using the Open Source Physics (OSP) Java library. The primary goal was to minimize overlaps and bounding box area (Packing) or ensure constrained placement within a predefined area (Floor-planning).


Task Implementation Summary
Task 1:
        The PackingSimulationExtended.java code simulates a system of circular particles of unit mass interacting under long-range attraction and short-range repulsion. Initial positions and sizes are randomly assigned, while velocities are initialized to zero. The Verlet algorithm is used to update positions and velocities iteratively. Air resistance is introduced to promote convergence by damping velocities. Repeated executions of this code with varied parameters revealed steady clustering behavior with minimal overlap when force magnitudes were properly balanced.


Tasks 2 and 3: 
        Particles were simulated under varying force parameters. Logging mechanisms recorded particle positions, bounding box dimensions, and system energy across time. These experiments showed that higher attraction magnitudes or lower repulsion strengths often resulted in increased overlap. When unit density was assumed instead of unit mass (implying mass scaled with size), heavier particles exhibited more inertia, reducing responsiveness to forces. The system thus required finer tuning to reach steady states, but once tuned, it produced more naturally spaced layouts.


Task 4:
        The FloorplanningSimulation.java code introduced a bounding area within which all particles must remain. The simulation logic ensured that particle positions did not escape the defined region, introducing additional boundary-handling conditions. Compared to the packing task, the final configurations were more compact and conformed clearly to the rectangular constraints. This made the simulations slower to converge, particularly with high particle densities, but the results were well-contained and maintained minimal overlap.


Tasks 5 and 6:
        To reduce the time complexity from , we introduced spatial partitioning in PackingSimulationOptimized.java and FloorplanningSimulationOptimized.java. By dividing the simulation space into square lattice cells, we localized force calculations. Interactions with distant cells were replaced with cumulative forces with particles that represented the center of mass. This optimization improved performance without significant loss in simulation accuracy. Comparisons of bounding box areas and particle layouts confirmed that optimized results closely matched unoptimized ones in both tasks, but with improved speed, especially for large N.
