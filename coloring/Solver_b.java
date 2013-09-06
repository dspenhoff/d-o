import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.Random;

public class Solver_b {
  
  int N, E;
  Graph g;
  int numUsedColors = 0;
  int numVColored = 0;
  int MAX_COLORS = 1000;
  int[] countColor = new int[MAX_COLORS];  
  Vertix[] vertices;    // do not assume *any* order (may be sorted or not)
  Color[] colors;       // do not assume *any* order (may be sorted or not)
  Random rng = new Random();
  Stopwatch timer;
  
  // parameters for simulated annealing
  boolean anneal = true;
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
      Solver_b solver = new Solver_b(args);
    } catch (IOException e) {
        e.printStackTrace();
    }
  }
    
  /**
  * Read the instance, solve it, and print the solution in the standard output
  */
  public Solver_b(String[] args) throws IOException {
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
        
    vertices = new Vertix[g.V()];
    for (int v = 0; v < g.V(); v++) {
      vertices[v] = new Vertix(v);
    }
    
    // greedy feasible solution
    
    /*
    NOTE: the quality of the final coloring varies depending on the sort order of vertices 
    and on which greedy method is used to start the process with a feasible coloring
    i used different combinations to get the 3 scores of 10; I believe the score of 7
    can be achieved on all cases regardless.
    */
    
    timer = new Stopwatch();

    Bag<Color> colorBag = new Bag<Color>();
    java.util.Arrays.sort(vertices);
    numUsedColors = 0;
    while (numVColored < g.V()) {
      Color kolor = new Color(numUsedColors++);
      colorBag.add(kolor);
      for (int v = 0; v < g.V(); v++) {
        if (vertices[v].color() != null) continue;  // v is already colored
        colorVertix(vertices[v], kolor);
      }
    }
    
    // put the colors in an array (for sorting, etc.)
    colors = new Color[colorBag.size()];
    int k = 0;
    for (Color color : colorBag) {
      colors[k++] = new Color(color);
    }
    
    System.out.println("elapsed time, color count = " + timer.elapsedTime() + ", " + colorCount());
    printColorCountViolations();    
    
/*  
    // alternate greedy
    java.util.Arrays.sort(vertices);
    Color[] kolors = new Color[1000];
    numUsedColors = 0;
    for (int v = 0; v < g.V(); v++) {
      boolean coloredVertix = false;
      for (int c = 0; c < numUsedColors; c++) {
        if (colorVertix(vertices[v], kolors[c])) {
          coloredVertix = true;
          break;
        } 
      }
      if (!coloredVertix) {
        kolors[numUsedColors] = new Color(numUsedColors);
        colorBag.add(kolors[numUsedColors]);
        numUsedColors++;
        colorVertix(vertices[v], kolors[numUsedColors-1]);          
      }
    }
*/ 
/*
    // yet another alternate greedy
    
    // open an initial set of colors (based on previous heuristic)
    colorBag = new Bag<Color>();
    numUsedColors *= 0.9;
    for (int i = 0; i < numUsedColors; i++) {
      colorBag.add(new Color(i));
    }
    
    uncolorAll();
    
    java.util.Arrays.sort(vertices);
    while (numVColored < g.V()) {
      for (Vertix v : vertices) {
        if (!v.isColored()) {
          boolean colored = false;
          for (Color c : colorBag) {
            if (colorVertix(v, c)) {
              colored = true;
              break;
            }
          }
          if (!colored) {
            // need a new color for v
            Color c = new Color(numUsedColors++);
            colorBag.add(c);
            colorVertix(v, c);          
          }
        }
      }
    }
    
    // put the colors in an array (for sorting, etc.)
    colors = new Color[colorBag.size()];
    k = 0;
    for (Color color : colorBag) {
      colors[k++] = new Color(color);
    }

    System.out.println("elapsed time, color count = " + timer.elapsedTime() + ", " + colorCount());
    printColorCountViolations();
    //printSolution();
 */   
        
    // improve greedy by eliminating a color?
    
    /*
    todo: try a kempe chain implementation
    */
    
    System.out.println("elapsed time, color count = " + timer.elapsedTime() + ", " + colorCount());
    printColorCountViolations();
    printSolution();
    //changeSingleColors();
    
    changeMultipleColors(100);
    
    System.out.println("elapsed time, color count = " + timer.elapsedTime() + ", " + colorCount());
    printColorCountViolations();
      
    numUsedColors = colorCount();
    //printColors();
    // display the solution
    printSolution();
    printViolations();

  }
  
  private void changeSingleColors() {
    // looks for opportunities to change the color of single vertices
    // in a way that tries to eliminate some used colors
    
    double Z = objective();
    System.out.println("objective = " + Z);
    
    java.util.Arrays.sort(vertices);
    int maxTries = Math.max(g.V(), 200);
    for (int tries = 0; tries < maxTries; tries++) {
      boolean solutionChanged = false;      
      java.util.Arrays.sort(colors);   // sorts colors by increasing count
      for (Color c : colors) {
        for (Vertix v : vertices) {
          if (c.count() == 0) break;     // skip "dead" colors
          if (v.color().color() == c.color()) {
            // v is a candidate to recolor
            // which non-adjacent color has highest count? 
            Color k = maxNeighborUnusedColor(v);
            if (k == null) continue;
            double deltac = - c.count() * c.count() + (c.count() - 1) * (c.count() - 1);
            double deltak = - k.count() * k.count() + (k.count() + 1) * (k.count() + 1);
            double delta = deltac + deltak;
            if (delta > 0 || anneal(delta)) {
              Z += deltac + deltak;
              colorVertix(v, k);          
              c.setUnused();
              if (!c.used()) break;
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
      Vertix v = vertices[rng.nextInt(g.V())];
      Vertix w = randomNeighbor(v);
      if (numNeighborsColor(v, w.color()) > 1) continue; 
      v.setChain();
      Bag<Vertix> chain = kempe(v.color(), w.color(), v);
      if (chain.size() == 0) continue;
      
      chain.add(v);
      
      // print the chain
      System.out.println("\nchain:");
      for (Vertix u : chain) {
        System.out.println("vertix, color = " + u.vertix() + ", " + u.color());
      }
      double delta = kempeDelta(chain, v.color(), w.color());
      System.out.println("delta = " + delta);
      
      if (delta > 0 || anneal(delta)) {
        // swap the chain colors
        //int sumv = 0;
        //int sumw = 0;
        for (Vertix u : chain) {
          u.clearChain();
          if (u.color().color() == v.color().color()) {
            // swap w for v
            colorVertix(u, w.color());
            colors[whichColor(v.color().color())].setUnused();
            //colors[v.color().color()].setUnused();
          } else {
            // swap v for w
            colorVertix(u, v.color());
            colors[whichColor(w.color().color())].setUnused();
            //colors[w.color().color()].setUnused();
          }
          u.clearChain();
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
  
  private Bag<Vertix> kempe(Color a, Color b, Vertix v) {
    // if possible, extends a kempe chain in colors a and b from vertix v
    // in the recursion, the chain objevt (bag) is instantiated when the stopping condition is reached
    // vertices are then added to the chain and the recursion unwind
    Bag<Vertix> chain = new Bag<Vertix>();
    for (int j : g.adj(v.vertix())) {
      Vertix w = vertices[whichVertix(j)];
      if (w.color().color() == b.color()) {
        if (w.chain()) continue;      // w is already in the chain
        if (!w.validColor(a)) break;
        w.setChain();
        chain = kempe(b, a, w);
        //w.clearChain();
        chain.add(w);
        break;        
      }
    }
    return chain;   // empty bag returned indicates stopping condition was reached
  }
  
  private double kempeDelta(Bag<Vertix> chain, Color a, Color b) {
    // return the delta to the objective if the colors in the kempe chain were swapped
    // the chain consists of vertices of the two given colors
    
    int suma = 0;
    int sumb = 0;
    for (Vertix v : chain) {
      if (v.color().color() == a.color()) suma++; 
      else                                sumb++;
    }
    
    int counta = a.count();
    int countb = b.count();
    int countaSwap = counta - suma + sumb;
    int countbSwap = countb - sumb + suma;
    double delta = - (counta * counta + countb * countb);
    return delta + (countaSwap * countaSwap + countbSwap * countbSwap);
  }
  
  
  private double objective() {
    double sum = 0;
    for (Color c : colors) sum += c.count() * c.count();
    return sum;
  }
    
  private Bag<Vertix> neighbors(Vertix v) {
    Bag<Vertix> neighborBag = new Bag<Vertix>();
    for (int u : g.adj(v.vertix())) {
      neighborBag.add(vertices[whichVertix(u)]);
    }
    return neighborBag;
  }
  
  private int numNeighborsColor(Vertix v, Color c) {
    // how many of vertix's neighbors are color?
    int count = 0;
    for (int i : g.adj(v.vertix())) {
      if (vertices[i].color().color() == c.color()) count++;
    }
    return count;
  }
  
  private Vertix randomNeighbor(Vertix v) {
    int n = rng.nextInt(v.degree());
    int j = 0;
    for (int i : g.adj(v.vertix())) if (j++ == n) return vertices[whichVertix(i)];
    return null;
  }
  
  private Color maxNeighborUnusedColor(Vertix v) {
    Bag<Vertix> neighborhood = neighbors(v);
    if (neighborhood.size() == numUsedColors-1) {
      return null;
    }
    boolean[] usedColors = new boolean[numUsedColors];
    for (Vertix w : neighborhood) {
      usedColors[w.color().color()] = true;
    }
    usedColors[v.color().color()] = true;
    
    int used = 0;
    for (int i = 0; i < numUsedColors; i++) if (usedColors[i]) used++;
    //System.out.println("used = " + used);
    
    int maxCount = 0;
    Color maxColor = null;
    for (int i = 0; i < numUsedColors; i++) {
      if (!usedColors[i]) {
        if (v.color().color() == i) continue;
        if (colors[whichColor(i)].count() > maxCount) {
          maxColor = colors[whichColor(i)];
          maxCount = maxColor.count();
        }
      }
    }
    //System.out.println("max color count = " + maxCount);
    return maxColor;
  }
  
  private void printColors() {
    for (Color c : colors) {
      System.out.println(c);
    }
  }

  private int whichVertix(int v) {
    // returns the index in vertices[] that represents the vertix "v"
    // method required when vertices[] has been sorted
    for (int i = 0; i < g.V(); i++) {
      if (vertices[i].vertix() == v) return i;
    }
    return -1;    // could be used for robust error checking
  }

  private int whichColor(int c) {
    // returns the index in colors[] that represents the color "c"
    // method required when colors[] has been sorted
    for (int i = 0; i < colors.length; i++) {
      if (colors[i].color() == c) return i;
    }
    return -1;    // could be used for robust error checking
  }
       
  private boolean colorVertix(Vertix v, Color c) {
    // colors vertix v with color c if possible (violates no constraint)
    // returns true for success, false for failure
    // note: will recolor if v is already colored
    for (Vertix w : neighbors(v)) if (w.color() == c) return false;  // not feasible
    if (v.color() == null) numVColored++;
    v.setColor(c);
    c.setUsed();
    return true;     
  }
  
  private int colorCount() {
    // how many colors were used?
    boolean[] usedColors = new boolean[numUsedColors];
    for (int v = 0; v < g.V(); v++) {
      usedColors[vertices[v].color().color()] = true;
    }
    int count = 0;
    for (int c = 0; c < colors.length; c++) {
      if (usedColors[c]) count++;
    }
    return count;    
  }
  
  private void printViolations() {
    boolean flag = false;
    String s = "Constraint violations: ";
    for (int v = 0; v < g.V(); v++) {
      for (int w : g.adj(vertices[v].vertix())) {
        for (int u = 0; u < g.V(); u++) {
          if (vertices[u].vertix() == w) {
            if (vertices[u].color() == vertices[v].color()) {
              flag = true;
              s += "(" + v + "," + u + ")";
              break;
            }
          } 
        }
      }
    }
    if (flag) System.out.println(s);
  }
  
  private void printSolution() {
    // print the solution in the specified output format        
    System.out.println(numUsedColors + " 0");
    String s = "";
    for (int v = 0; v < g.V(); v++) {
      for (int w = 0; w < g.V(); w++) {
        if (vertices[w].vertix() == v) {
          s += vertices[w].color().color() + " ";
          break;
        }
      }
    }
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
  
  
  private void printColorCountViolations() {
    // checks and prints inconsistencies betwen color counts and actual color assignments
    System.out.println("color count violations");
    int[] countingColors = new int[numUsedColors];
    for (int v = 0; v < g.V(); v++) {
      countingColors[vertices[v].color().color()]++;
    } 
    for (Color c : colors) {
      if (countingColors[c.color()] != c.count()) {
        System.out.println(c + " >> actual count = " + countingColors[c.color()]);
      }
    }  
  }
  
  private void uncolorAll() {
    for (Vertix v : vertices) v.setColor(null);
    numVColored = 0;
  }
  
  private class Vertix implements Comparable<Vertix> {
    private int vertix;
    private int degree;
    private Color color;
    private boolean chain;
    public Vertix(int v) {
      vertix = v;
      degree = 0;
      color = null;
      chain = false;
      for (int i : g.adj(vertix)) degree++;
    }
    public int compareTo(Vertix that) {
      // decreasing order
      if (this.degree > that.degree()) return -1;
      if (this.degree < that.degree()) return 1;
      return 0;
    }
    public int vertix() { return vertix; }
    public int degree() { return degree; }
    public void setColor(Color c) { color = c; }
    public Color color() { return color; }
    public boolean isColored() { return color != null; }
    public void setChain() { chain = true; }
    public void clearChain() { chain = false; }
    public boolean chain() { return chain; }
    public boolean validColor(Color c) {
      for (int i : g.adj(vertix)) {
        if (vertices[i].chain()) continue;
        if (vertices[i].color().color() == c.color()) return false;
      }
      return true;
    }
    public String toString() {
      String s = "vertix, degree, color: " + vertix + ", " + degree + ", " + color;
      s += "\nneighbors: ";
      for (int w : g.adj(vertix)) {
        s += w + " ";
      }
      return s;
    }
  }
  
  private class Color implements Comparable<Color> {
    private int color;
    private int count;
    private boolean used;
    public Color(int c) {
      color = c;
      count = 0;
      used = false;
    }
    public Color(Color c) {
      color = c.color();
      count = c.count();
      used = c.used();      
    }
    public int compareTo(Color that) {
      // increasing order
      if (this.count > that.count()) return 1;
      if (this.count < that.count()) return -1;
      return 0;      
    }
    public int color() { return color; }
    public void setColor(int c) { color = c; }
    public int count() { return count; }
    public void setCount(int n) { 
      count = n; 
      if (count == 0) used = false; 
    }
    public boolean used() { return used; }
    public void setUsed() {
      count++;
      used = true;
    }
    public void setUnused() {
      if (count > 0) count--;
      if (count == 0) used = false;
    }
    public String toString() {
      return "color, count = " + color + ", " + count;
    }
  }
  

}