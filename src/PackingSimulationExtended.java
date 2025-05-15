import java.awt.*;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import javax.swing.*;
import java.io.*;
import java.util.*;
import java.text.DecimalFormat;
import java.util.List;

public class PackingSimulationExtended extends JPanel {
    static class Particle {
        double x, y, vx = 0, vy = 0, ax = 0, ay = 0, radius;
        double mass; // Mass based on unit density
        Color color;
        boolean isOverlapping;

        Particle(double x, double y, double radius) {
            this.x = x;
            this.y = y;
            this.radius = radius;
            this.mass = Math.PI * radius * radius; // Area = mass (unit density)
            this.color = new Color(100, 100, 255, 150); // Semi-transparent
            this.isOverlapping = false;
        }
    }

    // Simulation parameters
    static double dt = 0.01;
    static int T = 15000; // Increased iterations
    static double damping = 0.97;
    static double k_repel = 15.0; // Stronger repulsion
    static double boundaryStiffness = 20.0;
    static List<Particle> particles = new ArrayList<>();
    static final int WIDTH = 1000, HEIGHT = 1000;
    static double minX, maxX, minY, maxY;
    static DecimalFormat df = new DecimalFormat("0.0000");
    static int overlapCount = 0;
    static String currentConfig = "default";

    public static void main(String[] args) throws IOException {
        new File("task3_logs").mkdirs();
        new File("task3_screenshots").mkdirs();

        // Parameter sets (dt, T, damping, k_repel, boundaryStiffness, particleCount)
        Object[][] configs = {
            {0.03, 12000, 0.96, 12.0, 15.0, 15},
            {0.02, 15000, 0.97, 15.0, 20.0, 15},
            {0.01, 20000, 0.98, 18.0, 25.0, 15}
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
            generateRandomParticles(particleCount);

            JFrame frame = new JFrame("Density Packing: " + currentConfig);
            PackingSimulationExtended panel = new PackingSimulationExtended();
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.add(panel);
            frame.pack();
            frame.setVisible(true);

            runSimulation(panel);
            frame.dispose();
        }
    }

    static void generateRandomParticles(int N) {
        Random rand = new Random();
        double spacing = 3.5; // Increased initial spacing
        int gridSize = (int)Math.ceil(Math.sqrt(N));
        
        for (int i = 0; i < N; i++) {
            int row = i / gridSize;
            int col = i % gridSize;
            double x = 4.0 + col * spacing;
            double y = 4.0 + row * spacing;
            double radius = 0.6 + rand.nextDouble() * 0.4; // 0.6-1.0 radius
            particles.add(new Particle(x, y, radius));
        }
        updateBoundaries();
    }

    static void runSimulation(PackingSimulationExtended panel) throws IOException {
        String logPrefix = "task3_logs/" + currentConfig;
        String screenshotPrefix = "task3_screenshots/" + currentConfig;

        try (BufferedWriter paramsLog = new BufferedWriter(new FileWriter(logPrefix + "_params.txt"));
             BufferedWriter energyLog = new BufferedWriter(new FileWriter(logPrefix + "_energy.csv"));
             BufferedWriter boundsLog = new BufferedWriter(new FileWriter(logPrefix + "_bounds.csv"))) {
       
            // Log parameters
            paramsLog.write(String.format(
                "Configuration: %s\n" +
                "dt: %.4f\nT: %d\ndamping: %.4f\nk_repel: %.4f\n" +
                "boundaryStiffness: %.4f\nParticles: %d\n\n",
                currentConfig, dt, T, damping, k_repel, boundaryStiffness, particles.size()));

            // Log initial state
            paramsLog.write("Initial positions (x y radius mass):\n");
            for (Particle p : particles) {
                paramsLog.write(String.format("%.4f %.4f %.4f %.4f\n", 
                    p.x, p.y, p.radius, p.mass));
            }

            // CSV headers
            energyLog.write("time,potential_energy,kinetic_energy,overlaps\n");
            boundsLog.write("time,width,height,area,density\n");

            // Main loop
            for (int t = 0; t < T; t++) {
                detectOverlaps();
                computeForces();
                updateParticles();
                updateBoundaries();

                if (t % 100 == 0 || t == T - 1) {
                    double time = t * dt;
                    double[] energies = calculateEnergies();
                    double totalArea = (maxX-minX)*(maxY-minY);
                    double totalMass = particles.stream().mapToDouble(p -> p.mass).sum();
                    double density = totalMass/totalArea;

                    // Log data
                    energyLog.write(String.format("%s,%s,%s,%d\n",
                        df.format(time), df.format(energies[0]), 
                        df.format(energies[1]), overlapCount));
                    boundsLog.write(String.format("%s,%s,%s,%s,%s\n",
                        df.format(time), df.format(maxX-minX),
                        df.format(maxY-minY), df.format(totalArea),
                        df.format(density)));

                    // Visual update
                    if (t % 500 == 0 || t == T - 1) {
                        panel.repaint();
                        try { Thread.sleep(10); } catch (InterruptedException e) {}
                    }
                }
            }

            // Save final state
            paramsLog.write("\nFinal positions (x y radius mass):\n");
            for (Particle p : particles) {
                paramsLog.write(String.format("%.4f %.4f %.4f %.4f\n", 
                    p.x, p.y, p.radius, p.mass));
            }
            
            double finalArea = (maxX-minX)*(maxY-minY);
            double finalDensity = particles.stream().mapToDouble(p -> p.mass).sum() / finalArea;
            paramsLog.write(String.format("\nFinal bounds: %.4f x %.4f\nArea: %.4f\nDensity: %.4f",
                maxX-minX, maxY-minY, finalArea, finalDensity));

            saveScreenshot(panel, screenshotPrefix + "_final.png");
        }
    }

    static void detectOverlaps() {
        overlapCount = 0;
        for (Particle p : particles) p.isOverlapping = false;

        for (int i = 0; i < particles.size(); i++) {
            Particle p1 = particles.get(i);
            for (int j = i+1; j < particles.size(); j++) {
                Particle p2 = particles.get(j);
                double dx = p2.x - p1.x;
                double dy = p2.y - p1.y;
                double dist = Math.sqrt(dx*dx + dy*dy);
                double minDist = p1.radius + p2.radius;

                if (dist < minDist) {
                    p1.isOverlapping = true;
                    p2.isOverlapping = true;
                    overlapCount++;
                }
            }
        }
    }

    static void computeForces() {
        // Reset accelerations
        for (Particle p : particles) {
            p.ax = 0;
            p.ay = 0;
        }

        // Particle interactions
        for (int i = 0; i < particles.size(); i++) {
            Particle p1 = particles.get(i);
            for (int j = i+1; j < particles.size(); j++) {
                Particle p2 = particles.get(j);
                double dx = p2.x - p1.x;
                double dy = p2.y - p1.y;
                double dist = Math.sqrt(dx*dx + dy*dy);
                double minDist = p1.radius + p2.radius;

                if (dist < minDist && dist > 1e-6) {
                    double overlap = minDist - dist;
                    double force = k_repel * Math.pow(overlap, 1.5) / p1.mass; // Mass-dependent
                    double fx = force * dx / dist;
                    double fy = force * dy / dist;

                    p1.ax -= fx;
                    p1.ay -= fy;
                    p2.ax += fx;
                    p2.ay += fy;
                }
            }
        }

        // Boundary forces
        for (Particle p : particles) {
            double leftOverlap = minX - (p.x - p.radius);
            if (leftOverlap > 0) p.ax += boundaryStiffness * leftOverlap / p.mass;
            
            double rightOverlap = (p.x + p.radius) - maxX;
            if (rightOverlap > 0) p.ax -= boundaryStiffness * rightOverlap / p.mass;
            
            double bottomOverlap = minY - (p.y - p.radius);
            if (bottomOverlap > 0) p.ay += boundaryStiffness * bottomOverlap / p.mass;
            
            double topOverlap = (p.y + p.radius) - maxY;
            if (topOverlap > 0) p.ay -= boundaryStiffness * topOverlap / p.mass;
        }
    }

    static void updateParticles() {
        for (Particle p : particles) {
            // Verlet integration with position correction
            double prevX = p.x;
            double prevY = p.y;
            
            p.x += p.vx * dt + 0.5 * p.ax * dt * dt;
            p.y += p.vy * dt + 0.5 * p.ay * dt * dt;
            
            // Stronger damping for overlapping particles
            double localDamping = p.isOverlapping ? damping * 0.8 : damping;
            p.vx = (p.vx + 0.5 * p.ax * dt) * localDamping;
            p.vy = (p.vy + 0.5 * p.ay * dt) * localDamping;
            
            // Position correction
            if (p.isOverlapping) {
                p.x = 0.8 * prevX + 0.2 * p.x;
                p.y = 0.8 * prevY + 0.2 * p.y;
                p.vx *= 0.7;
                p.vy *= 0.7;
            }
        }
    }

    static double[] calculateEnergies() {
        double potential = 0;
        double kinetic = 0;
        
        // Particle interactions
        for (int i = 0; i < particles.size(); i++) {
            Particle p1 = particles.get(i);
            kinetic += 0.5 * p1.mass * (p1.vx*p1.vx + p1.vy*p1.vy);
            
            for (int j = i+1; j < particles.size(); j++) {
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
        
        // Boundary interactions
        for (Particle p : particles) {
            double leftOverlap = minX - (p.x - p.radius);
            if (leftOverlap > 0) potential += 0.5 * boundaryStiffness * leftOverlap * leftOverlap;
            
            double rightOverlap = (p.x + p.radius) - maxX;
            if (rightOverlap > 0) potential += 0.5 * boundaryStiffness * rightOverlap * rightOverlap;
            
            double bottomOverlap = minY - (p.y - p.radius);
            if (bottomOverlap > 0) potential += 0.5 * boundaryStiffness * bottomOverlap * bottomOverlap;
            
            double topOverlap = (p.y + p.radius) - maxY;
            if (topOverlap > 0) potential += 0.5 * boundaryStiffness * topOverlap * topOverlap;
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
        
        // Draw particles as empty circles with mass-dependent size
        g2d.setStroke(new BasicStroke(1.5f));
        for (Particle p : particles) {
            int x = (int)(p.x * scale + offsetX);
            int y = (int)(p.y * scale + offsetY);
            int r = (int)(p.radius * scale);
            
            // Draw particle
            g2d.setColor(p.color);
            g2d.drawOval(x - r, y - r, 2 * r, 2 * r);
            
            // Indicate overlaps
            if (p.isOverlapping) {
                g2d.setColor(new Color(255, 0, 0, 100));
                g2d.fillOval(x - r, y - r, 2 * r, 2 * r);
            }
            
            // Show mass text for larger particles
            if (p.radius > 0.8) {
                g2d.setColor(Color.BLACK);
                g2d.setFont(new Font("Arial", Font.PLAIN, 10));
                g2d.drawString(String.format("%.2f", p.mass), x - 10, y + 4);
            }
        }
        
        // Draw bounding box
        g2d.setColor(Color.RED);
        g2d.setStroke(new BasicStroke(2));
        g2d.drawRect(
            (int)(minX * scale + offsetX),
            (int)(minY * scale + offsetY),
            (int)((maxX - minX) * scale),
            (int)((maxY - minY) * scale)
        );
        
        // Display info
        double totalMass = particles.stream().mapToDouble(p -> p.mass).sum();
        double area = (maxX-minX)*(maxY-minY);
        double density = totalMass/area;
        
        g2d.setColor(Color.BLACK);
        g2d.setFont(new Font("Arial", Font.BOLD, 14));
        g2d.drawString("Config: " + currentConfig, 10, 20);
        g2d.drawString(String.format("Particles: %d  Overlaps: %d", 
            particles.size(), overlapCount), 10, 40);
        g2d.drawString(String.format("Density: %.4f", density), 10, 60);
        g2d.drawString(String.format("Area: %.2f", area), 10, 80);
    }

    @Override
    public Dimension getPreferredSize() {
        return new Dimension(WIDTH, HEIGHT);
    }
}