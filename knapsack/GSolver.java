import java.io.*;
import java.util.List;
import java.util.ArrayList;

public class GSolver {
  
  private int capacity;
  private int N;
  private Item[] items;
      
  /**
  * The main class
  */
  public static void main(String[] args) {
    try {
      GSolver solver = new GSolver(args);
    } catch (IOException e) {
        e.printStackTrace();
    }
  }
  
  
  /**
   * greedy solver for the knapsack problem
   */
  public GSolver(String[] args) throws IOException {
    
    // get the temp file name
    String fileName = null;
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
    capacity = Integer.parseInt(firstLine[1]);
    items = new Item[N];
    
    for(int i = 0; i < N; i++){
      String line = lines.get(i+1);
      String[] parts = line.split("\\s+");

      int v = Integer.parseInt(parts[0]);
      int w = Integer.parseInt(parts[1]);
      items[i] = new Item(v, w, i);
    }

    // a greedy algorithm that sorts the items by "bang for the buck"
    // and takes items that fit in bang-order

    // sort the items based on "bang"
    
    java.util.Arrays.sort(items);

    // 1st pass - basic greedy heuristic
    int sumW = 0;
    int sumV = 0;
    int lastSelected = -1;
    for (int i = 0; i < N; i++) {
      Item item = items[i];
      if ((sumW + item.weight()) < capacity) {
        item.setSelected();
        sumW += item.weight();
        sumV += item.value();
        lastSelected = i;
      }
    }
    
    System.out.println(value());
        
    
  }
  
  private void doRandomSwaps(int n, int K) {
    // attempts K times to randomly swap out an item in current solution for best greedy replacement(s)
    
    for (int k = 0; k < K; k++) {
 
      trySwap(i, j);
    }
  }
  
  private void randomlSwap(int n, int K) {
    // attempts K times to randomly swap out (up to) n items in current solution for best greedy replacement(s)
    int z = value();
    
    for int k = 0; k < K; k++) {
      
    }
    
  }
  
  private boolean trySwap()
  
  private int value() {
    int z = 0;
    for (int i = 0; i < N; i++) if (items[i].isSelected()) z += items[i].value();
    return z;
  }
  
  
  /*
  *  data structure for the items
  */
  
  private class Item implements Comparable<Item> {
    int value;
    int weight;
    double bang;
    int index;
    boolean selected;
    public Item( int v, int w, int i) {
      value = v;
      weight = w;
      index = i;
      selected = false;
      bang = (double) value / weight;        
    }
    
    public int value() { return value; }
    public int weight() { return weight; }
    public double bang() { return bang; }
    public void setSelected() { selected = true; }
    public void setNotSelected() { selected = false; }
    public int index() { return index; }
    public boolean isSelected() { return selected; }
    
    public int compareTo(Item that) {
      // note: returns reverse order compare
      if (this.bang() > that.bang()) {
        return -1;
      } else if (this.bang() < that.bang()) {
        return 1;
      }
      return 0;
    }
    
    public String toString() {
      String s = "";
      s += "Item: " + index;
      s += ", value: " + value;
      s += ", weight: " + weight;
      s += ", bang: " + bang;
      return s;
    }
  }
  
  
}