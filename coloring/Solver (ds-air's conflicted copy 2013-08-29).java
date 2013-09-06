import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.BitSet;

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
  Vertix[] vertices;    // do not assume *any* order (may be sorted or not)
  Color[] colors;       // do not assume *any* order (may be sorted or not)
    
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
        
    vertices = new Vertix[g.V()];
    for (int v = 0; v < g.V(); v++) {
      vertices[v] = new Vertix(v);
    }
    
    Bag<Color> colorBag = new Bag<Color>();

    // greedy feasible solution
    
    /*
    NOTE: the quality of the final coloring varies depending on the sort order of vertices 
    and on which greedy method is used to start the process with a feasible coloring
    i used different combinations to get the 3 scores of 10; I believe the score of 7
    can be achieved on all cases regardless.
    */
/*
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
*/    
    // alternate greedy
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
    
    // improve greedy by eliminating a color?
    
    colors = new Color[colorBag.size()];
    int k = 0;
    for (Color color : colorBag) {
      colors[k++] = new Color(color);
    }
    
    /*
    todo: try a kempe chain implementation
    */
      
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
            Color newColor = maxNeighborUnusedColor(v);
            if (newColor == null) continue;
            Color oldColor = v.color();
            colorVertix(v, newColor);          
            oldColor.setUnused();
            if (!oldColor.used()) break;
            solutionChanged = true;
          }
        }
      }
      if (!solutionChanged) break;
    }

    // how many colors were used?
    boolean[] usedColors = new boolean[numUsedColors];
    for (int v = 0; v < g.V(); v++) {
      usedColors[vertices[v].color().color()] = true;
    }
    int colorCount = 0;
    for (int c = 0; c < colors.length; c++) {
      if (usedColors[c]) colorCount++;
    }
    numUsedColors = colorCount;
    
    // display the solution
    printSolution();
    printViolations();

  }
    
  private Bag<Vertix> neighbors(Vertix v) {
    Bag<Vertix> neighborBag = new Bag<Vertix>();
    for (int u : g.adj(v.vertix())) {
      neighborBag.add(vertices[whichVertix(u)]);
    }
    return neighborBag;
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
  
  private class Vertix implements Comparable<Vertix> {
    private int vertix;
    private int degree;
    Color color;
    public Vertix(int v) {
      vertix = v;
      degree = 0;
      for (int i : g.adj(v)) degree++;
    }
    public int compareTo(Vertix that) {
      // increasing order
      if (this.degree > that.degree()) return -1;
      if (this.degree < that.degree()) return 1;
      return 0;
    }
    public int vertix() { return vertix; }
    public int degree() { return degree; }
    public void setColor(Color c) { color = c; }
    public Color color() { return color; }
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
    public boolean used() { return used; }
    public void setUsed() {
      count++;
      used = true;
    }
    public void setUnused() {
      count--;
      if (count == 0) used = false;
    }
    public String toString() {
      return "color, count = " + color + ", " + count;
    }
  }
  

}