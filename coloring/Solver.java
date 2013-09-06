import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.Random;

/**
 * The class <code>Solver</code> is an implementation of a greedy algorithm to solve the knapsack problem.
 *
*/

public class Solver {
  
  int N, E;
  Graph g;
  int numUsedColors = 0;
  int numVColored = 0;
  int MAX_COLORS = 1000;
  int[] countColor = new int[MAX_COLORS];  
  int[] vertices;     // color of vertix j; -1 means no color assigned
  int[] degree;       // degree of vertix j
  boolean[] inChain;  // vertix j is in the kempe chain
  int[] colors;       // the number of times color i is used; 0 means not used
  Random rng = new Random();
  Stopwatch timer;
  
  // parameters for simulated annealing
  boolean anneal = false;
  double temperatureStart = 1000;
  double temperature = temperatureStart;
  double temperatureMin = 0.001;
  double temperatureDelta = 1 - temperatureMin;
  double temperatureCount = 0;
  double temperatureChange = 1000;
    
  /**
  * The main class
  */
  public static void main(String[] args) {
    try {
      Solver solver = new Solver(args);
    } catch (IOException e) {
        e.printStackTrace();
    }
  }
    
  /**
  * Read the instance, solve it, and print the solution in the standard output
  */
  public Solver(String[] args) throws IOException {
    String fileName = null;
        
    // get the temp file name
    for(String arg : args){
      if(arg.startsWith("-file=")){
        fileName = arg.substring(6);
      } 
    }
    if(fileName == null)
      return;
        
    // read the lines out of the file
    List<String> lines = new ArrayList<String>();

    BufferedReader input =  new BufferedReader(new FileReader(fileName));
    try {
      String line = null;
      while (( line = input.readLine()) != null){
        lines.add(line);
      }
    }
    finally {
      input.close();
    }
              
    // parse the data in the file
    String[] firstLine = lines.get(0).split("\\s+");
    N = Integer.parseInt(firstLine[0]);
    E = Integer.parseInt(firstLine[1]);
    g = new Graph(N);

    for(int i = 0; i < E; i++){
      String line = lines.get(i+1);
      String[] parts = line.split("\\s+");
      int v = Integer.parseInt(parts[0]);
      int w = Integer.parseInt(parts[1]);
      g.addEdge(v, w);
    }  
        
    vertices = new int[g.V()];
    degree = new int[g.V()];
    inChain = new boolean[g.V()];
    for (int j = 0; j < g.V(); j++) { 
      vertices[j] = -1;
      for (int jj : g.adj(j)) degree[j]++;
    }
    
    // greedy feasible solution
    
    
    timer = new Stopwatch();

    colors = new int[g.V()];    // max colors = number of vertices

    numUsedColors = 0;
    while (numVColored < g.V()) {
      int kolor = numUsedColors++;
      for (int j = 0; j < g.V(); j++) {
        if (vertices[j] < 0) colorVertix(j, kolor);   // try to color uncolored vertices
      }
    }
    
    System.out.println("elapsed time = " + timer.elapsedTime());
    System.out.println("greedy assignment:");
    printSolution();
    printViolations();
    
    changeMultipleColors(1000000);
    
    System.out.println("elapsed time = " + timer.elapsedTime());
    System.out.println("chain swaps:");
    printSolution();
    printViolations();

  }
  
  private void changeSingleColors() {
    // looks for opportunities to change the color of single vertices
    // in a way that tries to eliminate some used colors
    
    double Z = objective();
    System.out.println("objective = " + Z);
    

    int maxTries = Math.max(g.V(), 200);
    for (int tries = 0; tries < maxTries; tries++) {
      boolean solutionChanged = false;      

      for (int c = 0; c < numUsedColors; c++) {
        for (int v : vertices) {
          if (colors[c] == 0) break;     // skip "dead" colors
          if (vertices[v] == c) {
            // v is a candidate to recolor
            // which non-adjacent color has highest count? 
            int k = maxNeighborUnusedColor(v);
            if (k == -1) continue;
            double deltac = - colors[c] * colors[c] + (colors[c] - 1) * (colors[c] - 1);
            double deltak = - colors[k] * colors[k] + (colors[k] + 1) * (colors[k] + 1);
            double delta = deltac + deltak;
            if (delta > 0 || anneal(delta)) {
              Z += deltac + deltak;
              colorVertix(v, k);          
              colors[c]--;
              if (colors[c] == 0) break;
              solutionChanged = true;              
            }
          }
        }
      }
      if (!solutionChanged) break;
    }
    System.out.println("objective = " + objective());
  }
  
  private void changeMultipleColors(int n) {
    // looks for opportunintes to change colors in a chain of vertices
    // eg, kempe chain
    // todo
    
    // randomly select a vertix to change
    // randomly select an adjacent vertix
    // build the associated kempe chain
    // check the delta and swap the chain
    
    double Z = objective();
    System.out.println("objective = " + Z);
    
    for (int j = 0; j < n; j++) {
      int v = vertices[rng.nextInt(g.V())];
      int w = randomNeighbor(v);
      if (numNeighborsColor(v, colors[w]) > 1) continue; 
      inChain[v] = true;
      Bag<Integer> chain = kempe(colors[v], colors[w], v);
      chain.add(v);
      if (chain.size() == 1) continue;
      
      
      // print the chain
      System.out.println("\nchain:");
      for (int u : chain) {
        System.out.println("vertix, color = " + u + ", " + vertices[u]);
      }
      double delta = kempeDelta(chain, vertices[v], vertices[w]);
      System.out.println("delta = " + delta);
      
      if (delta > 0 || anneal(delta)) {
        // swap the chain colors
        //int sumv = 0;
        //int sumw = 0;
        System.out.println(">>> swapping chain");
        for (int u : chain) {
          inChain[u] = false;
          int cv = colors[v];
          int cw = colors[w];
          if (colors[u] == cv) {
            // swap w for v
            colors[cv]--;
            vertices[u] = cw;
            colors[cw]++;
            System.out.println("swapping w for v: " + u + " from " + colors[u] + " to " + colors[w]);
            System.out.println(colors[w]);
            System.out.println(colors[v]);
            //colors[v.color().color()].setUnused();
          } else {
            // swap v for w
            colors[cw]--;
            vertices[u] = cv;
            colors[cv]++;
            System.out.println("swapping v for w: " + u + " from " + colors[u] + " to " + colors[v]);
            System.out.println(colors[v]);
            System.out.println(colors[w]);            
            //colors[w.color().color()].setUnused();
          }
          inChain[u] = false;
        }
        
        //Color cv = colors[whichColor(v.color().color())];
        //Color cw = colors[whichColor(w.color().color())];
        //cv.setCount(cv.count() - sumv + sumw);
        //cw.setCount(cw.count() - sumv + sumw);
        Z += delta;
        //System.out.println("Z = " + Z);
        
      }
    }
    System.out.println("Z = " + Z);
    System.out.println("objective = " + objective());
    
  }
  
  private Bag<Integer> kempe(int colora, int colorb, int vertix) {
    // if possible, extends a kempe chain in colors a and b from vertix v
    // in the recursion, the chain objevt (bag) is instantiated when the stopping condition is reached
    // vertices are then added to the chain and the recursion unwind
    Bag<Integer> chain = new Bag<Integer>();
    for (int j : g.adj(vertix)) {
      if (vertices[j] == colorb) {
        if (inChain[j]) continue;      // w is already in the chain
        if (validColor(j, colora)) {
          inChain[j] = true;
          chain = kempe(colorb, colora, j);
          chain.add(j);
          break;          
        } 
      }
    }
    return chain;   // empty bag returned indicates stopping condition was reached
  }
  
  private double kempeDelta(Bag<Integer> chain, int colora, int colorb) {
    // return the delta to the objective if the colors in the kempe chain were swapped
    // the chain consists of vertices of the two given colors
    
    int suma = 0;
    int sumb = 0;
    for (int v : chain) {
      if (vertices[v] == colora) suma++; 
      else                       sumb++;
    }
    
    int counta = colors[colora];
    int countb = colors[colorb];
    
    int countaSwap = counta - suma + sumb;
    int countbSwap = countb - sumb + suma;
    double delta = - (counta * counta + countb * countb);
    delta += (countaSwap * countaSwap + countbSwap * countbSwap);
    System.out.println("suma, sumb, counta, countb, delta = " + suma + " " + sumb + " " + counta + " " + countb + " " + delta);
    return delta;
  }
  
  
  private double objective() {
    double sum = 0;
    for (int i = 0; i < numUsedColors; i++) sum += colors[i] * colors[i];
    return sum;
  }
  
  private Bag<Integer> neighbors(int vertix) {
    Bag<Integer> neighborBag = new Bag<Integer>();
    for (int u : g.adj(vertix)) neighborBag.add(u);
    return neighborBag;
  }  
  
  private int numNeighborsColor(int vertix, int color) {
    // how many of vertix's neighbors are assigned color?
    int count = 0;
    for (int j : g.adj(vertix)) {
      if (vertices[j] == color) count++;
    }
    return count;
  }
  
  private int randomNeighbor(int vertix) {
    int n = rng.nextInt(degree[vertix]);
    int j = 0;
    for (int i : g.adj(vertix)) if (j++ == n) return i;
    return -1;
  }
  
  private int maxNeighborUnusedColor(int vertix) {
    Bag<Integer> neighbs = neighbors(vertix);
    if (neighbs.size() == numUsedColors-1) return -1;
    boolean[] usedColors = new boolean[numUsedColors];
    for (int j : neighbs) usedColors[vertices[j]] = true;
    usedColors[vertices[vertix]] = true;
    
    int used = 0;
    for (int i = 0; i < numUsedColors; i++) if (usedColors[i]) used++;
    
    int maxCount = 0;
    int maxColor = -1;
    for (int i = 0; i < numUsedColors; i++) {
      if (!usedColors[i]) {
        if (vertices[vertix] == i) continue;
        if (colors[i] > maxCount) {
          maxColor = colors[i];
          maxCount = maxColor;
        }
      }
    }

    return maxColor;
  }
  
  private void printColors() {
    for (int i = 0; i < numUsedColors; i++) System.out.println("color, count = " + i + ", " + colors[i]);
  }

  private boolean colorVertix(int vertix, int color) {
    // colors vertix with color if possible (violates no constraint)
    // returns true for success, false for failure
    // note: will recolor if vertix is already colored
    for (int j : neighbors(vertix)) if (vertices[j] == color) return false;  // not feasible
    if (vertices[vertix] == -1) numVColored++;
    vertices[vertix] = color;
    colors[color]++;
    return true;     
  }
  
  private int colorCount() {
    // how many colors were used?
    boolean[] usedColors = new boolean[numUsedColors];
    for (int j = 0; j < g.V(); j++) usedColors[vertices[j]] = true;
    int count = 0;
    for (int i = 0; i < numUsedColors; i++) if (usedColors[i]) count++;
    return count;    
  }
  
  private void printViolations() {
    boolean flag = false;
    String s = "Constraint violations: ";
    for (int v = 0; v < g.V(); v++) {
      for (int w : g.adj(v)) {
        if (w == v) continue;
        if (vertices[w] == vertices[v]) {
          flag = true;
          s += "(" + v + "," + w + ")";
          break;
        }
      }
    }
    if (flag) System.out.println(s);
  }
  
  private void printSolution() {
    // print the solution in the specified output format        
    System.out.println(numUsedColors + " 0");
    String s = "";
    for (int j = 0; j < g.V(); j++) s += vertices[j] + " ";
    System.out.println(s);
  }
  
  private boolean anneal(double delta) {
    if (!anneal || temperature < temperatureMin) return false;
    boolean accept = (rng.nextDouble() < Math.exp(-delta / temperature));
    if (accept && temperatureCount++ > temperatureChange) {
      temperature *= temperatureDelta;
      temperatureCount = 0;
    }
    return accept;
  }
  
  public boolean validColor(int vertix, int color) {
    for (int i : g.adj(vertix)) {
      if (inChain[i]) continue;
      if (vertices[i] == color) return false;
    }
    return true;
  }
  
  private void printColorCountViolations() {
    // checks and prints inconsistencies betwen color counts and actual color assignments
    System.out.println("color count violations");
    int[] countingColors = new int[numUsedColors];
    for (int v = 0; v < g.V(); v++) countingColors[vertices[v]]++;
    for (int i = 0; i < numUsedColors; i++) {
      if (countingColors[i] != colors[i]) {
        System.out.println(i + ", " + colors[i] + " >> actual count = " + countingColors[i]);
      }
    }  
  }
  
  private void uncolorAll() {
    for (int v : vertices) vertices[v] = -1;
    numVColored = 0;
  }

}