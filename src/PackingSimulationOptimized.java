import java.awt.*;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import javax.swing.*;
import java.io.*;
import java.util.*;
import java.text.DecimalFormat;
import java.util.List;

public class PackingSimulationOptimized extends JPanel {
    static class Particle {
        double x, y, vx = 0, vy = 0, ax = 0, ay = 0, radius;
        double mass;
        Color color;
        boolean isOverlapping;

        Particle(double x, double y, double radius) {
            this.x = x;
            this.y = y;
            this.radius = radius;
            this.mass = Math.PI * radius * radius;
            this.color = new Color(100, 100, 255, 150);
            this.isOverlapping = false;
        }
    }

    static class Cell {
        List<Particle> particles = new ArrayList<>();
        double centerX, mass;
        double cx, cy; // Center coordinates (geometric or center of mass)

        Cell(double cx, double cy) {
            this.cx = cx;
            this.cy = cy;
        }

        void computeCenter(boolean useCenterOfMass) {
            if (particles.isEmpty()) {
                mass = 0;
                return;
            }
            if (useCenterOfMass) {
                double totalMass = 0;
                double sumX = 0, sumY = 0;
                for (Particle p : particles) {
                    sumX += p.x * p.mass;
                    sumY += p.y * p.mass;
                    totalMass += p.mass;
                }
                this.mass = totalMass;
                this.cx = sumX / totalMass;
                this.cy = sumY / totalMass;
            } else {
                this.mass = particles.stream().mapToDouble(p -> p.mass).sum();
                // cx, cy already set to geometric center
            }
        }
    }

    // Simulation parameters
    static double dt = 0.01;
    static int T = 15000;
    static double damping = 0.97;
    static double k_repel = 15.0;
    static double boundaryStiffness = 20.0;
    static List<Particle> particles = new ArrayList<>();
    static final int WIDTH = 1000, HEIGHT = 1000;
    static double minX, maxX, minY, maxY;
    static DecimalFormat df = new DecimalFormat("0.0000");
    static int overlapCount = 0;
    static String currentConfig = "default";
    static Cell[][] grid;
    static double cellSize;

    public static void main(String[] args) throws IOException {
        new File("task5_logs").mkdirs();
        new File("task5_screenshots").mkdirs();

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
            loadInitialConditions("task3_logs/config_" + (i+1) + "_params.txt", particleCount);
            initializeGrid();

            JFrame frame = new JFrame("Optimized Packing: " + currentConfig);
            PackingSimulationOptimized panel = new PackingSimulationOptimized();
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.add(panel);
            frame.pack();
            frame.setVisible(true);

            runSimulation(panel);
            frame.dispose();
        }
    }

    static void loadInitialConditions(String filename, int particleCount) throws IOException {
        particles.clear();
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line;
            boolean readingPositions = false;
            int count = 0;
            while ((line = br.readLine()) != null && count < particleCount) {
                if (line.startsWith("Initial positions")) {
                    readingPositions = true;
                    continue;
                }
                if (readingPositions && !line.trim().isEmpty()) {
                    String[] parts = line.trim().split("\\s+");
                    if (parts.length >= 3) {
                        double x = Double.parseDouble(parts[0]);
                        double y = Double.parseDouble(parts[1]);
                        double radius = Double.parseDouble(parts[2]);
                        particles.add(new Particle(x, y, radius));
                        count++;
                    }
                }
            }
        }
        updateBoundaries();
    }

    static void initializeGrid() {
        double range = Math.max(maxX - minX, maxY - minY);
        cellSize = 2.0; // Adjust based on max radius
        int gridSize = (int)Math.ceil(range / cellSize);
        grid = new Cell[gridSize][gridSize];
        for (int i = 0; i < gridSize; i++) {
            for (int j = 0; j < gridSize; j++) {
                double cx = minX + (i + 0.5) * cellSize;
                double cy = minY + (j + 0.5) * cellSize;
                grid[i][j] = new Cell(cx, cy);
            }
        }
        assignParticlesToCells();
    }

    static void assignParticlesToCells() {
        for (Cell[] row : grid) {
            for (Cell cell : row) {
                cell.particles.clear();
            }
        }
        for (Particle p : particles) {
            int i = (int)((p.x - minX) / cellSize);
            int j = (int)((p.y - minY) / cellSize);
            if (i >= 0 && i < grid.length && j >= 0 && j < grid[0].length) {
                grid[i][j].particles.add(p);
            }
        }
        for (Cell[] row : grid) {
            for (Cell cell : row) {
                cell.computeCenter(false); // Geometric center initially
            }
        }
    }

    static void runSimulation(PackingSimulationOptimized panel) throws IOException {
        String logPrefix = "task5_logs/" + currentConfig;
        String screenshotPrefix = "task5_screenshots/" + currentConfig;

        try (BufferedWriter paramsLog = new BufferedWriter(new FileWriter(logPrefix + "_params.txt"));
             BufferedWriter energyLog = new BufferedWriter(new FileWriter(logPrefix + "_energy.csv"));
             BufferedWriter boundsLog = new BufferedWriter(new FileWriter(logPrefix + "_bounds.csv"))) {

            paramsLog.write(String.format(
                "Configuration: %s\n" +
                "dt: %.4f\nT: %d\ndamping: %.4f\nk_repel: %.4f\n" +
                "boundaryStiffness: %.4f\nParticles: %d\nCellSize: %.4f\n\n",
                currentConfig, dt, T, damping, k_repel, boundaryStiffness, particles.size(), cellSize));

            paramsLog.write("Initial positions (x y radius mass):\n");
            for (Particle p : particles) {
                paramsLog.write(String.format("%.4f %.4f %.4f %.4f\n", 
                    p.x, p.y, p.radius, p.mass));
            }

            energyLog.write("time,potential_energy,kinetic_energy,overlaps\n");
            boundsLog.write("time,width,height,area,density\n");

            for (int t = 0; t < T; t++) {
                detectOverlaps();
                computeForces();
                updateParticles();
                assignParticlesToCells();
                updateBoundaries();

                if (t % 100 == 0 || t == T - 1) {
                    double time = t * dt;
                    double[] energies = calculateEnergies();
                    double totalArea = (maxX-minX)*(maxY-minY);
                    double totalMass = particles.stream().mapToDouble(p -> p.mass).sum();
                    double density = totalMass/totalArea;

                    energyLog.write(String.format("%s,%s,%s,%d\n",
                        df.format(time), df.format(energies[0]), 
                        df.format(energies[1]), overlapCount));
                    boundsLog.write(String.format("%s,%s,%s,%s,%s\n",
                        df.format(time), df.format(maxX-minX),
                        df.format(maxY-minY), df.format(totalArea),
                        df.format(density)));

                    if (t % 500 == 0 || t == T - 1) {
                        panel.repaint();
                        try { Thread.sleep(10); } catch (InterruptedException e) {}
                    }
                }
            }

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

        for (int i = 0; i < grid.length; i++) {
            for (int j = 0; j < grid[0].length; j++) {
                Cell cell = grid[i][j];
                for (int di = -1; di <= 1; di++) {
                    for (int dj = -1; dj <= 1; dj++) {
                        int ni = i + di;
                        int nj = j + dj;
                        if (ni >= 0 && ni < grid.length && nj >= 0 && nj < grid[0].length) {
                            Cell neighbor = grid[ni][nj];
                            for (Particle p1 : cell.particles) {
                                for (Particle p2 : neighbor.particles) {
                                    if (p1 != p2) {
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
                        }
                    }
                }
            }
        }
    }

    static void computeForces() {
        for (Particle p : particles) {
            p.ax = 0;
            p.ay = 0;
        }

        for (int i = 0; i < grid.length; i++) {
            for (int j = 0; j < grid[0].length; j++) {
                Cell cell = grid[i][j];
                cell.computeCenter(true); // Center of mass for interactions
                for (Particle p1 : cell.particles) {
                    // Interactions within same cell and adjacent cells
                    for (int di = -1; di <= 1; di++) {
                        for (int dj = -1; dj <= 1; dj++) {
                            int ni = i + di;
                            int nj = j + dj;
                            if (ni >= 0 && ni < grid.length && nj >= 0 && nj < grid[0].length) {
                                Cell neighbor = grid[ni][nj];
                                if (di == 0 && dj == 0) {
                                    // Same cell: individual interactions
                                    for (Particle p2 : neighbor.particles) {
                                        if (p1 != p2) {
                                            computeParticleInteraction(p1, p2);
                                        }
                                    }
                                } else {
                                    // Adjacent cell: super-particle interaction
                                    if (neighbor.mass > 0) {
                                        computeSuperParticleInteraction(p1, neighbor);
                                    }
                                }
                            }
                        }
                    }
                    // Non-adjacent cells: super-particle at geometric center
                    for (int ni = 0; ni < grid.length; ni++) {
                        for (int nj = 0; nj < grid[0].length; nj++) {
                            if (Math.abs(ni - i) > 1 || Math.abs(nj - j) > 1) {
                                Cell neighbor = grid[ni][nj];
                                neighbor.computeCenter(false); // Geometric center
                                if (neighbor.mass > 0) {
                                    computeSuperParticleInteraction(p1, neighbor);
                                }
                            }
                        }
                    }
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

    static void computeParticleInteraction(Particle p1, Particle p2) {
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
            p2.ax += fx;
            p2.ay += fy;
        }
    }

    static void computeSuperParticleInteraction(Particle p1, Cell cell) {
        double dx = cell.cx - p1.x;
        double dy = cell.cy - p1.y;
        double dist = Math.sqrt(dx*dx + dy*dy);
        double minDist = p1.radius; // Super-particle has no radius

        if (dist < minDist && dist > 1e-6) {
            double overlap = minDist - dist;
            double force = k_repel * Math.pow(overlap, 1.5) / p1.mass;
            double fx = force * dx / dist;
            double fy = force * dy / dist;

            p1.ax -= fx;
            p1.ay -= fy;
        }
    }

    static void updateParticles() {
        for (Particle p : particles) {
            double prevX = p.x;
            double prevY = p.y;
            
            p.x += p.vx * dt + 0.5 * p.ax * dt * dt;
            p.y += p.vy * dt + 0.5 * p.ay * dt * dt;
            
            double localDamping = p.isOverlapping ? damping * 0.8 : damping;
            p.vx = (p.vx + 0.5 * p.ax * dt) * localDamping;
            p.vy = (p.vy + 0.5 * p.ay * dt) * localDamping;
            
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
        
        for (int i = 0; i < grid.length; i++) {
            for (int j = 0; j < grid[0].length; j++) {
                Cell cell = grid[i][j];
                cell.computeCenter(true);
                for (Particle p1 : cell.particles) {
                    kinetic += 0.5 * p1.mass * (p1.vx*p1.vx + p1.vy*p1.vy);
                    for (int di = -1; di <= 1; di++) {
                        for (int dj = -1; dj <= 1; dj++) {
                            int ni = i + di;
                            int nj = j + dj;
                            if (ni >= 0 && ni < grid.length && nj >= 0 && nj < grid[0].length) {
                                Cell neighbor = grid[ni][nj];
                                if (di == 0 && dj == 0) {
                                    for (Particle p2 : neighbor.particles) {
                                        if (p1 != p2) {
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
                                } else if (neighbor.mass > 0) {
                                    double dx = neighbor.cx - p1.x;
                                    double dy = neighbor.cy - p1.y;
                                    double dist = Math.sqrt(dx*dx + dy*dy);
                                    double minDist = p1.radius;
                                    if (dist < minDist) {
                                        double overlap = minDist - dist;
                                        potential += 0.5 * k_repel * Math.pow(overlap, 2.5);
                                    }
                                }
                            }
                        }
                    }
                    for (int ni = 0; ni < grid.length; ni++) {
                        for (int nj = 0; nj < grid[0].length; nj++) {
                            if (Math.abs(ni - i) > 1 || Math.abs(nj - j) > 1) {
                                Cell neighbor = grid[ni][nj];
                                neighbor.computeCenter(false);
                                if (neighbor.mass > 0) {
                                    double dx = neighbor.cx - p1.x;
                                    double dy = neighbor.cy - p1.y;
                                    double dist = Math.sqrt(dx*dx + dy*dy);
                                    double minDist = p1.radius;
                                    if (dist < minDist) {
                                        double overlap = minDist - dist;
                                        potential += 0.5 * k_repel * Math.pow(overlap, 2.5);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

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
        
        g2d.setColor(Color.WHITE);
        g2d.fillRect(0, 0, getWidth(), getHeight());
        
        double scaleX = (getWidth() - 40) / (maxX - minX);
        double scaleY = (getHeight() - 40) / (maxY - minY);
        double scale = Math.min(scaleX, scaleY);
        double offsetX = (getWidth() - (maxX-minX)*scale)/2 - minX*scale;
        double offsetY = (getHeight() - (maxY-minY)*scale)/2 - minY*scale;
        
        // Draw grid
        g2d.setColor(new Color(200, 200, 200));
        g2d.setStroke(new BasicStroke(0.5f));
        for (int i = 0; i <= grid.length; i++) {
            int x = (int)(minX * scale + offsetX + i * cellSize * scale);
            g2d.drawLine(x, (int)(minY * scale + offsetY), x, (int)(maxY * scale + offsetY));
        }
        for (int j = 0; j <= grid[0].length; j++) {
            int y = (int)(minY * scale + offsetY + j * cellSize * scale);
            g2d.drawLine((int)(minX * scale + offsetX), y, (int)(maxX * scale + offsetX), y);
        }

        // Draw particles
        g2d.setStroke(new BasicStroke(1.5f));
        for (Particle p : particles) {
            int x = (int)(p.x * scale + offsetX);
            int y = (int)(p.y * scale + offsetY);
            int r = (int)(p.radius * scale);
            
            g2d.setColor(p.color);
            g2d.drawOval(x - r, y - r, 2 * r, 2 * r);
            
            if (p.isOverlapping) {
                g2d.setColor(new Color(255, 0, 0, 100));
                g2d.fillOval(x - r, y - r, 2 * r, 2 * r);
            }
            
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