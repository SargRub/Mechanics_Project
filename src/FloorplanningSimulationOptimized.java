import java.awt.*;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import javax.swing.*;
import java.io.*;
import java.util.*;
import java.text.DecimalFormat;
import java.util.List;

public class FloorplanningSimulationOptimized extends JPanel {
    static class Particle {
        double x, y, vx = 0, vy = 0, ax = 0, ay = 0, radius;
        double mass;
        Color color;
        boolean isOverlapping;
        boolean isFixed;

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

    static class Cell {
        List<Particle> particles = new ArrayList<>();
        double mass;
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
    static double boundaryStiffness = 25.0;
    static double targetArea = 100.0;
    static List<Particle> particles = new ArrayList<>();
    static final int WIDTH = 1000, HEIGHT = 1000;
    static double minX, maxX, minY, maxY;
    static DecimalFormat df = new DecimalFormat("0.0000");
    static int overlapCount = 0;
    static String currentConfig = "default";
    static Cell[][] grid;
    static double cellSize;

    public static void main(String[] args) throws IOException {
        new File("task6_logs").mkdirs();
        new File("task6_screenshots").mkdirs();

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
            int movableParticleCount = (int)configs[i][5];

            particles.clear();
            loadInitialConditions("task4_logs/config_" + (i+1) + "_params.txt", movableParticleCount);
            initializeGrid();

            JFrame frame = new JFrame("Optimized Floorplanning: " + currentConfig);
            FloorplanningSimulationOptimized panel = new FloorplanningSimulationOptimized();
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.add(panel);
            frame.pack();
            frame.setVisible(true);

            runSimulation(panel);
            frame.dispose();
        }
    }

    static void loadInitialConditions(String filename, int movableParticleCount) throws IOException {
    particles.clear();
    int expectedParticleCount = movableParticleCount + 4; // 4 fixed particles
    try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
        String line;
        boolean readingPositions = false;
        int count = 0;
        while ((line = br.readLine()) != null && count < expectedParticleCount) {
            if (line.startsWith("Initial positions")) {
                readingPositions = true;
                continue;
            }
            if (readingPositions && !line.trim().isEmpty()) {
                // Check if line is a header or non-data (e.g., "Final positions")
                if (line.trim().startsWith("Final") || line.trim().startsWith("Target")) {
                    break; // Stop reading at next section
                }
                String[] parts = line.trim().split("\\s+");
                if (parts.length >= 5) {
                    try {
                        double x = Double.parseDouble(parts[0]);
                        double y = Double.parseDouble(parts[1]);
                        double radius = Double.parseDouble(parts[2]);
                        boolean isFixed = Boolean.parseBoolean(parts[4]);
                        particles.add(new Particle(x, y, radius, isFixed));
                        count++;
                    } catch (NumberFormatException e) {
                        // Skip malformed lines
                        continue;
                    }
                }
            }
        }
        if (count < expectedParticleCount) {
            System.err.println("Warning: Loaded " + count + " particles, expected " + expectedParticleCount);
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

    static void runSimulation(FloorplanningSimulationOptimized panel) throws IOException {
        String logPrefix = "task6_logs/" + currentConfig;
        String screenshotPrefix = "task6_screenshots/" + currentConfig;

        try (BufferedWriter paramsLog = new BufferedWriter(new FileWriter(logPrefix + "_params.txt"));
             BufferedWriter energyLog = new BufferedWriter(new FileWriter(logPrefix + "_energy.csv"));
             BufferedWriter boundsLog = new BufferedWriter(new FileWriter(logPrefix + "_bounds.csv"))) {

            paramsLog.write(String.format(
                "Configuration: %s\nTarget Area: %.4f\n" +
                "dt: %.4f\nT: %d\ndamping: %.4f\nk_repel: %.4f\n" +
                "boundaryStiffness: %.4f\nMovable Particles: %d\nCellSize: %.4f\n\n",
                currentConfig, targetArea, dt, T, damping, 
                k_repel, boundaryStiffness, particles.stream().filter(p -> !p.isFixed).count(), cellSize));

            paramsLog.write("Initial positions (x y radius mass isFixed):\n");
            for (Particle p : particles) {
                paramsLog.write(String.format("%.4f %.4f %.4f %.4f %b\n", 
                    p.x, p.y, p.radius, p.mass, p.isFixed));
            }

            energyLog.write("time,potential_energy,kinetic_energy,overlaps\n");
            boundsLog.write("time,width,height,area,actual_area,density\n");

            for (int t = 0; t < T; t++) {
                detectOverlaps();
                computeForces();
                updateParticles();
                assignParticlesToCells();
                updateBoundaries();

                if (t % 100 == 0 || t == T - 1) {
                    double time = t * dt;
                    double[] energies = calculateEnergies();
                    double actualArea = (maxX-minX)*(maxY-minY);
                    double totalMass = particles.stream().filter(p -> !p.isFixed)
                                              .mapToDouble(p -> p.mass).sum();
                    double density = totalMass/targetArea;

                    energyLog.write(String.format("%s,%s,%s,%d\n",
                        df.format(time), df.format(energies[0]), 
                        df.format(energies[1]), overlapCount));
                    boundsLog.write(String.format("%s,%s,%s,%s,%s,%s\n",
                        df.format(time), df.format(maxX-minX),
                        df.format(maxY-minY), df.format(targetArea),
                        df.format(actualArea), df.format(density)));

                    if (t % 500 == 0 || t == T - 1) {
                        panel.repaint();
                        try { Thread.sleep(10); } catch (InterruptedException e) {}
                    }
                }
            }

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
                                if (p1.isFixed) continue;
                                for (Particle p2 : neighbor.particles) {
                                    if (p1 != p2) {
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
                        }
                    }
                }
            }
        }
    }

    static void computeForces() {
        for (Particle p : particles) {
            if (p.isFixed) continue;
            p.ax = 0;
            p.ay = 0;
        }

        for (int i = 0; i < grid.length; i++) {
            for (int j = 0; j < grid[0].length; j++) {
                Cell cell = grid[i][j];
                cell.computeCenter(true); // Center of mass for interactions
                for (Particle p1 : cell.particles) {
                    if (p1.isFixed) continue;
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
            
            if (p.x - p.radius < left) {
                double overlap = left - (p.x - p.radius);
                p.ax += boundaryStiffness * overlap / p.mass;
            }
            if (p.x + p.radius > right) {
                double overlap = (p.x + p.radius) - right;
                p.ax -= boundaryStiffness * overlap / p.mass;
            }
            if (p.y - p.radius < bottom) {
                double overlap = bottom - (p.y - p.radius);
                p.ay += boundaryStiffness * overlap / p.mass;
            }
            if (p.y + p.radius > top) {
                double overlap = (p.y + p.radius) - top;
                p.ay -= boundaryStiffness * overlap / p.mass;
            }
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
            if (!p2.isFixed) {
                p2.ax += fx;
                p2.ay += fy;
            }
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
            if (p.isFixed) continue;
            
            double prevX = p.x;
            double prevY = p.y;
            
            p.x += p.vx * dt + 0.5 * p.ax * dt * dt;
            p.y += p.vy * dt + 0.5 * p.ay * dt * dt;
            
            double localDamping = p.isOverlapping ? damping * 0.7 : damping;
            p.vx = (p.vx + 0.5 * p.ax * dt) * localDamping;
            p.vy = (p.vy + 0.5 * p.ay * dt) * localDamping;
            
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
        
        for (int i = 0; i < grid.length; i++) {
            for (int j = 0; j < grid[0].length; j++) {
                Cell cell = grid[i][j];
                cell.computeCenter(true);
                for (Particle p1 : cell.particles) {
                    if (p1.isFixed) continue;
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

        double areaSide = Math.sqrt(targetArea);
        double centerX = 10.0;
        double centerY = 10.0;
        double left = centerX - areaSide/2;
        double right = centerX + areaSide/2;
        double bottom = centerY - areaSide/2;
        double top = centerY + areaSide/2;

        for (Particle p : particles) {
            if (p.isFixed) continue;
            
            if (p.x - p.radius < left) {
                double overlap = left - (p.x - p.radius);
                potential += 0.5 * boundaryStiffness * overlap * overlap;
            }
            if (p.x + p.radius > right) {
                double overlap = (p.x + p.radius) - right;
                potential += 0.5 * boundaryStiffness * overlap * overlap;
            }
            if (p.y - p.radius < bottom) {
                double overlap = bottom - (p.y - p.radius);
                potential += 0.5 * boundaryStiffness * overlap * overlap;
            }
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

        // Draw target area
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
            
            g2d.setColor(p.color);
            g2d.drawOval(x - r, y - r, 2 * r, 2 * r);
            
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