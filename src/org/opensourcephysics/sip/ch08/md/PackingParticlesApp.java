package org.opensourcephysics.sip.ch08.md;

import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;

public class PackingParticlesApp extends AbstractSimulation {
    PackingParticles particles = new PackingParticles();
    DisplayFrame display = new DisplayFrame("x", "y", "Particle Packing");

    public void initialize() {
        particles.N = control.getInt("Number of particles");
        particles.Lx = control.getDouble("System width Lx");
        particles.Ly = control.getDouble("System height Ly");
        particles.dt = control.getDouble("Time step dt");
        particles.repulsionStrength = control.getDouble("Repulsion strength");
        particles.damping = control.getDouble("Damping coefficient");
        
        particles.initialize();
        display.addDrawable(particles);
        display.setPreferredMinMax(0, particles.Lx, 0, particles.Ly);
        display.setSquareAspect(true);
    }

    public void doStep() {
        particles.step();
        display.setMessage("Step: " + particles.steps);
        display.repaint();
    }

    public void reset() {
        control.setValue("Number of particles", 20);
        control.setAdjustableValue("System width Lx", 20.0);
        control.setAdjustableValue("System height Ly", 20.0);
        control.setAdjustableValue("Time step dt", 0.05);
        control.setValue("Repulsion strength", 50.0);
        control.setValue("Damping coefficient", 0.1);
        display.setSquareAspect(true);
    }

    public static void main(String[] args) {
        SimulationControl.createApp(new PackingParticlesApp());
    }
}