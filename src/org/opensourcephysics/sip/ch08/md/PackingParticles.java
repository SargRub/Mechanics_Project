package org.opensourcephysics.sip.ch08.md;

import java.awt.*;
import java.util.Arrays;
import org.opensourcephysics.display.*;
import org.opensourcephysics.numerics.*;

public class PackingParticles implements Drawable, ODE {
    // Simulation state
    public double[] state;  // [x1, vx1, y1, vy1, x2, vx2, ... t]
    public double[] ax, ay; // Accelerations
    public double[] radius; // Particle radii
    
    // Simulation parameters
    public int N = 20;      // Number of particles
    public double Lx = 20;  // System width
    public double Ly = 20;  // System height
    public double dt = 0.05; // Time step
    
    // Physics parameters
    public double repulsionStrength = 50.0;
    public double damping = 0.1;
    
    // Visualization parameters
    public double minRadius = 0.4;
    public double maxRadius = 0.8;
    public Color particleColor = Color.BLUE;
    
    // Simulation control
    public int steps = 0;
    public double t = 0;
    Verlet odeSolver = new Verlet(this);

    public void initialize() {
        // Initialize arrays
        state = new double[1 + 4*N]; // +1 for time
        ax = new double[N];
        ay = new double[N];
        radius = new double[N];
        
        // Set random positions and sizes
        setRandomPositionsAndSizes();
        
        // Initialize velocities to zero
        Arrays.fill(state, 0);
        
        // Initialize solver
        odeSolver.setStepSize(dt);
        steps = 0;
        t = 0;
    }

    public void setRandomPositionsAndSizes() {
        double padding = maxRadius * 2;
        for(int i = 0; i < N; i++) {
            radius[i] = minRadius + (maxRadius - minRadius) * Math.random();
            state[4*i] = padding + (Lx - 2*padding) * Math.random(); // x
            state[4*i+2] = padding + (Ly - 2*padding) * Math.random(); // y
        }
    }

    public void computeAcceleration() {
        Arrays.fill(ax, 0);
        Arrays.fill(ay, 0);
        
        for(int i = 0; i < N; i++) {
            for(int j = i+1; j < N; j++) {
                double dx = state[4*i] - state[4*j];
                double dy = state[4*i+2] - state[4*j+2];
                double r = Math.sqrt(dx*dx + dy*dy);
                double minDist = radius[i] + radius[j];
                
                if(r < minDist) {
                    double f = repulsionStrength * (minDist - r)/r;
                    ax[i] += f * dx;
                    ay[i] += f * dy;
                    ax[j] -= f * dx;
                    ay[j] -= f * dy;
                }
            }
            
            // Add damping
            ax[i] -= damping * state[4*i+1];
            ay[i] -= damping * state[4*i+3];
        }
    }

    public void step() {
        odeSolver.step();
        enforceBoundaries();
        t += dt;
        steps++;
    }

    private void enforceBoundaries() {
        for(int i = 0; i < N; i++) {
            double r = radius[i];
            state[4*i] = Math.max(r, Math.min(Lx - r, state[4*i]));
            state[4*i+2] = Math.max(r, Math.min(Ly - r, state[4*i+2]));
        }
    }

    public void draw(DrawingPanel panel, Graphics g) {
        if(state == null) return;
        
        g.setColor(particleColor);
        int pxRadius, pyRadius;
        
        for(int i = 0; i < N; i++) {
            pxRadius = Math.abs(panel.xToPix(radius[i]) - panel.xToPix(0));
            pyRadius = Math.abs(panel.yToPix(radius[i]) - panel.yToPix(0));
            
            int xpix = panel.xToPix(state[4*i]) - pxRadius;
            int ypix = panel.yToPix(state[4*i+2]) - pyRadius;
            
            g.fillOval(xpix, ypix, 2*pxRadius, 2*pyRadius);
        }
    }

    // ODE interface methods
    public void getRate(double[] state, double[] rate) {
        if(odeSolver.getRateCounter() == 1) {
            computeAcceleration();
        }
        
        for(int i = 0; i < N; i++) {
            rate[4*i] = state[4*i+1];   // dx/dt = vx
            rate[4*i+2] = state[4*i+3]; // dy/dt = vy
            rate[4*i+1] = ax[i];        // dvx/dt = ax
            rate[4*i+3] = ay[i];        // dvy/dt = ay
        }
        rate[4*N] = 1; // dt/dt = 1
    }

    public double[] getState() {
        return state;
    }
}