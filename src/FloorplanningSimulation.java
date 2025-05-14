import java.awt.*;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import javax.swing.*;
import java.io.*;
import java.util.*;
import java.text.DecimalFormat;
import java.util.List;

public class FloorplanningSimulation extends JPanel {
    static class Particle {
        double x, y, vx = 0, vy = 0, ax = 0, ay = 0, radius;
        double mass; // For density calculations
        Color color;
        boolean isOverlapping;
        boolean isFixed; // For fixed boundary particles

        Particle(double x, double y, double radius, boolean fixed) {
            this.x = x;
            this.y = y;
            this.radius = radius;
            this.mass = Math.PI * radius * radius;
            this.color = fixed ? new Color(200, 50, 50) : new Color(50, 50, 200);
            this.isOverlapping = false;
            this.isFixed = fixed;
        }
    }

    // Simulation parameters
    static double dt = 0.01;
    static int T = 15000;
    static double damping = 0.97;
    static double k_repel = 15.0;
    static double boundaryStiffness = 25.0;
    static double targetArea = 100.0; // Target area A for floorplanning
    static List<Particle> particles = new ArrayList<>();
    static final int WIDTH = 1000, HEIGHT = 1000;
    static double minX, maxX, minY, maxY;
    static DecimalFormat df = new DecimalFormat("0.0000");
    static int overlapCount = 0;
    static String currentConfig = "default";

    public static void main(String[] args) throws IOException {
        new File("task4_logs").mkdirs();
        new File("task4_screenshots").mkdirs();

        // Parameter sets (dt, T, damping, k_repel, boundaryStiffness, particleCount)
        Object[][] configs = {
            {0.03, 12000, 0.96, 12.0, 20.0, 12},
            {0.02, 15000, 0.97, 15.0, 25.0, 12},
            {0.01, 20000, 0.98, 18.0, 30.0, 12}
        };

        for (int i = 0; i < configs.length; i++) {
            currentConfig = "config_" + (i+1);
            dt = (double)configs[i][0];
            T = (int)configs[i][1];
            damping = (double)configs[i][2];
            k_repel = (double)configs[i][3];
            boundaryStiffness = (double)configs[i][4];
            int particleCount = (int)configs[i][5];

            particles.clear();
            initializeFloorplan(particleCount);

            JFrame frame = new JFrame("Floorplanning: " + currentConfig);
            FloorplanningSimulation panel = new FloorplanningSimulation();
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.add(panel);
            frame.pack();
            frame.setVisible(true);

            runSimulation(panel);
            frame.dispose();
        }
    }

    static void initializeFloorplan(int movableParticles) {
        // Create fixed boundary particles (rectangle area)
        double areaSide = Math.sqrt(targetArea);
        double centerX = 10.0;
        double centerY = 10.0;
        
        // Fixed boundary particles (4 corners)
        particles.add(new Particle(centerX - areaSide/2, centerY - areaSide/2, 0.5, true));
        particles.add(new Particle(centerX + areaSide/2, centerY - areaSide/2, 0.5, true));
        particles.add(new Particle(centerX - areaSide/2, centerY + areaSide/2, 0.5, true));
        particles.add(new Particle(centerX + areaSide/2, centerY + areaSide/2, 0.5, true));

        // Add movable particles inside the area
        Random rand = new Random();
        for (int i = 0; i < movableParticles; i++) {
            double x = centerX + (rand.nextDouble() - 0.5) * areaSide * 0.8;
            double y = centerY + (rand.nextDouble() - 0.5) * areaSide * 0.8;
            double radius = 0.5 + rand.nextDouble() * 0.5;
            particles.add(new Particle(x, y, radius, false));
        }
        updateBoundaries();
    }

    static void runSimulation(FloorplanningSimulation panel) throws IOException {
        String logPrefix = "task4_logs/" + currentConfig;
        String screenshotPrefix = "task4_screenshots/" + currentConfig;

        try (BufferedWriter paramsLog = new BufferedWriter(new FileWriter(logPrefix + "_params.txt"));
             BufferedWriter energyLog = new BufferedWriter(new FileWriter(logPrefix + "_energy.csv"));
             BufferedWriter boundsLog = new BufferedWriter(new FileWriter(logPrefix + "_bounds.csv"))) {

            // Log parameters
            paramsLog.write(String.format(
                "Configuration: %s\nTarget Area: %.4f\n" +
                "dt: %.4f\nT: %d\ndamping: %.4f\nk_repel: %.4f\n" +
                "boundaryStiffness: %.4f\nMovable Particles: %d\n\n",
                currentConfig, targetArea, dt, T, damping, 
                k_repel, boundaryStiffness, particles.stream().filter(p -> !p.isFixed).count()));

            // Log initial state
            paramsLog.write("Initial positions (x y radius mass isFixed):\n");
            for (Particle p : particles) {
                paramsLog.write(String.format("%.4f %.4f %.4f %.4f %b\n", 
                    p.x, p.y, p.radius, p.mass, p.isFixed));
            }

            // CSV headers
            energyLog.write("time,potential_energy,kinetic_energy,overlaps\n");
            boundsLog.write("time,width,height,area,actual_area,density\n");

            // Main loop
            for (int t = 0; t < T; t++) {
                detectOverlaps();
                computeForces();
                updateParticles();
                updateBoundaries();

                if (t % 100 == 0 || t == T - 1) {
                    double time = t * dt;
                    double[] energies = calculateEnergies();
                    double actualArea = (maxX-minX)*(maxY-minY);
                    double totalMass = particles.stream().filter(p -> !p.isFixed)
                                              .mapToDouble(p -> p.mass).sum();
                    double density = totalMass/targetArea;

                    // Log data
                    energyLog.write(String.format("%s,%s,%s,%d\n",
                        df.format(time), df.format(energies[0]), 
                        df.format(energies[1]), overlapCount));
                    boundsLog.write(String.format("%s,%s,%s,%s,%s,%s\n",
                        df.format(time), df.format(maxX-minX),
                        df.format(maxY-minY), df.format(targetArea),
                        df.format(actualArea), df.format(density)));

                    // Visual update
                    if (t % 500 == 0 || t == T - 1) {
                        panel.repaint();
                        try { Thread.sleep(10); } catch (InterruptedException e) {}
                    }
                }
            }

            // Save final state
            paramsLog.write("\nFinal positions (x y radius mass isFixed):\n");
            for (Particle p : particles) {
                paramsLog.write(String.format("%.4f %.4f %.4f %.4f %b\n", 
                    p.x, p.y, p.radius, p.mass, p.isFixed));
            }
            
            double finalActualArea = (maxX-minX)*(maxY-minY);
            double finalDensity = particles.stream().filter(p -> !p.isFixed)
                                         .mapToDouble(p -> p.mass).sum() / targetArea;
            paramsLog.write(String.format(
                "\nTarget Area: %.4f\nActual Area: %.4f\nDensity: %.4f",
                targetArea, finalActualArea, finalDensity));

            saveScreenshot(panel, screenshotPrefix + "_final.png");
        }
    }

    static void detectOverlaps() {
        overlapCount = 0;
        for (Particle p : particles) p.isOverlapping = false;

        for (int i = 0; i < particles.size(); i++) {
            Particle p1 = particles.get(i);
            if (p1.isFixed) continue;
            
            for (int j = 0; j < particles.size(); j++) {
                if (i == j) continue;
                Particle p2 = particles.get(j);
                double dx = p2.x - p1.x;
                double dy = p2.y - p1.y;
                double dist = Math.sqrt(dx*dx + dy*dy);
                double minDist = p1.radius + p2.radius;

                if (dist < minDist) {
                    p1.isOverlapping = true;
                    overlapCount++;
                    if (!p2.isFixed) p2.isOverlapping = true;
                }
            }
        }
    }

    static void computeForces() {
        // Reset accelerations
        for (Particle p : particles) {
            if (p.isFixed) continue;
            p.ax = 0;
            p.ay = 0;
        }

        // Particle interactions
        for (int i = 0; i < particles.size(); i++) {
            Particle p1 = particles.get(i);
            if (p1.isFixed) continue;
            
            for (int j = 0; j < particles.size(); j++) {
                if (i == j) continue;
                Particle p2 = particles.get(j);
                double dx = p2.x - p1.x;
                double dy = p2.y - p1.y;
                double dist = Math.sqrt(dx*dx + dy*dy);
                double minDist = p1.radius + p2.radius;

                if (dist < minDist && dist > 1e-6) {
                    double overlap = minDist - dist;
                    double force = k_repel * Math.pow(overlap, 1.5) / p1.mass;
                    double fx = force * dx / dist;
                    double fy = force * dy / dist;

                    p1.ax -= fx;
                    p1.ay -= fy;
                    if (!p2.isFixed) {
                        p2.ax += fx;
                        p2.ay += fy;
                    }
                }
            }
        }

        // Boundary constraints (stronger for floorplanning)
        double areaSide = Math.sqrt(targetArea);
        double centerX = 10.0;
        double centerY = 10.0;
        double left = centerX - areaSide/2;
        double right = centerX + areaSide/2;
        double bottom = centerY - areaSide/2;
        double top = centerY + areaSide/2;

        for (Particle p : particles) {
            if (p.isFixed) continue;
            
            // Left boundary
            if (p.x - p.radius < left) {
                double overlap = left - (p.x - p.radius);
                p.ax += boundaryStiffness * overlap / p.mass;
            }
            // Right boundary
            if (p.x + p.radius > right) {
                double overlap = (p.x + p.radius) - right;
                p.ax -= boundaryStiffness * overlap / p.mass;
            }
            // Bottom boundary
            if (p.y - p.radius < bottom) {
                double overlap = bottom - (p.y - p.radius);
                p.ay += boundaryStiffness * overlap / p.mass;
            }
            // Top boundary
            if (p.y + p.radius > top) {
                double overlap = (p.y + p.radius) - top;
                p.ay -= boundaryStiffness * overlap / p.mass;
            }
        }
    }

    static void updateParticles() {
        for (Particle p : particles) {
            if (p.isFixed) continue;
            
            // Verlet integration with position correction
            double prevX = p.x;
            double prevY = p.y;
            
            p.x += p.vx * dt + 0.5 * p.ax * dt * dt;
            p.y += p.vy * dt + 0.5 * p.ay * dt * dt;
            
            // Stronger damping for overlapping particles
            double localDamping = p.isOverlapping ? damping * 0.7 : damping;
            p.vx = (p.vx + 0.5 * p.ax * dt) * localDamping;
            p.vy = (p.vy + 0.5 * p.ay * dt) * localDamping;
            
            // Position correction
            if (p.isOverlapping) {
                p.x = 0.7 * prevX + 0.3 * p.x;
                p.y = 0.7 * prevY + 0.3 * p.y;
                p.vx *= 0.5;
                p.vy *= 0.5;
            }
        }
    }

    static double[] calculateEnergies() {
        double potential = 0;
        double kinetic = 0;
        
        // Particle interactions
        for (int i = 0; i < particles.size(); i++) {
            Particle p1 = particles.get(i);
            if (p1.isFixed) continue;
            
            kinetic += 0.5 * p1.mass * (p1.vx*p1.vx + p1.vy*p1.vy);
            
            for (int j = 0; j < particles.size(); j++) {
                if (i == j) continue;
                Particle p2 = particles.get(j);
                double dx = p2.x - p1.x;
                double dy = p2.y - p1.y;
                double dist = Math.sqrt(dx*dx + dy*dy);
                double minDist = p1.radius + p2.radius;
                
                if (dist < minDist) {
                    double overlap = minDist - dist;
                    potential += 0.5 * k_repel * Math.pow(overlap, 2.5);
                }
            }
        }
        
        // Boundary constraints
        double areaSide = Math.sqrt(targetArea);
        double centerX = 10.0;
        double centerY = 10.0;
        double left = centerX - areaSide/2;
        double right = centerX + areaSide/2;
        double bottom = centerY - areaSide/2;
        double top = centerY + areaSide/2;

        for (Particle p : particles) {
            if (p.isFixed) continue;
            
            // Left boundary
            if (p.x - p.radius < left) {
                double overlap = left - (p.x - p.radius);
                potential += 0.5 * boundaryStiffness * overlap * overlap;
            }
            // Right boundary
            if (p.x + p.radius > right) {
                double overlap = (p.x + p.radius) - right;
                potential += 0.5 * boundaryStiffness * overlap * overlap;
            }
            // Bottom boundary
            if (p.y - p.radius < bottom) {
                double overlap = bottom - (p.y - p.radius);
                potential += 0.5 * boundaryStiffness * overlap * overlap;
            }
            // Top boundary
            if (p.y + p.radius > top) {
                double overlap = (p.y + p.radius) - top;
                potential += 0.5 * boundaryStiffness * overlap * overlap;
            }
        }
        
        return new double[]{potential, kinetic};
    }

    static void updateBoundaries() {
        minX = minY = Double.MAX_VALUE;
        maxX = maxY = -Double.MAX_VALUE;
        for (Particle p : particles) {
            minX = Math.min(minX, p.x - p.radius);
            maxX = Math.max(maxX, p.x + p.radius);
            minY = Math.min(minY, p.y - p.radius);
            maxY = Math.max(maxY, p.y + p.radius);
        }
    }

    static void saveScreenshot(JPanel panel, String filename) throws IOException {
        BufferedImage img = new BufferedImage(
            panel.getWidth(), panel.getHeight(), BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = img.createGraphics();
        panel.paint(g2d);
        ImageIO.write(img, "png", new File(filename));
        g2d.dispose();
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2d = (Graphics2D)g;
        
        // White background
        g2d.setColor(Color.WHITE);
        g2d.fillRect(0, 0, getWidth(), getHeight());
        
        // Calculate scaling
        double scaleX = (getWidth() - 40) / (maxX - minX);
        double scaleY = (getHeight() - 40) / (maxY - minY);
        double scale = Math.min(scaleX, scaleY);
        double offsetX = (getWidth() - (maxX-minX)*scale)/2 - minX*scale;
        double offsetY = (getHeight() - (maxY-minY)*scale)/2 - minY*scale;
        
        // Draw target area (rectangle)
        double areaSide = Math.sqrt(targetArea);
        g2d.setColor(new Color(255, 200, 200));
        g2d.fillRect(
            (int)((10.0 - areaSide/2) * scale + offsetX),
            (int)((10.0 - areaSide/2) * scale + offsetY),
            (int)(areaSide * scale),
            (int)(areaSide * scale)
        );
        
        // Draw particles
        g2d.setStroke(new BasicStroke(1.5f));
        for (Particle p : particles) {
            int x = (int)(p.x * scale + offsetX);
            int y = (int)(p.y * scale + offsetY);
            int r = (int)(p.radius * scale);
            
            // Draw particle
            g2d.setColor(p.color);
            g2d.drawOval(x - r, y - r, 2 * r, 2 * r);
            
            // Indicate overlaps
            if (p.isOverlapping && !p.isFixed) {
                g2d.setColor(new Color(255, 0, 0, 100));
                g2d.fillOval(x - r, y - r, 2 * r, 2 * r);
            }
        }
        
        // Draw actual bounding box
        g2d.setColor(Color.RED);
        g2d.setStroke(new BasicStroke(2));
        g2d.drawRect(
            (int)(minX * scale + offsetX),
            (int)(minY * scale + offsetY),
            (int)((maxX - minX) * scale),
            (int)((maxY - minY) * scale)
        );
        
        // Display info
        double actualArea = (maxX-minX)*(maxY-minY);
        double movableMass = particles.stream().filter(p -> !p.isFixed)
                                    .mapToDouble(p -> p.mass).sum();
        double density = movableMass/targetArea;
        
        g2d.setColor(Color.BLACK);
        g2d.setFont(new Font("Arial", Font.BOLD, 14));
        g2d.drawString("Config: " + currentConfig, 10, 20);
        g2d.drawString(String.format("Target Area: %.2f  Actual: %.2f", 
            targetArea, actualArea), 10, 40);
        g2d.drawString(String.format("Density: %.4f", density), 10, 60);
        g2d.drawString(String.format("Overlaps: %d", overlapCount), 10, 80);
    }

    @Override
    public Dimension getPreferredSize() {
        return new Dimension(WIDTH, HEIGHT);
    }
}