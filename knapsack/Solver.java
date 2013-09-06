import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.BitSet;

/**
 * The class <code>Solver</code> is an implementation of a greedy algorithm to solve the knapsack problem.
 *
 */
public class Solver {
  
  int capacity;
  int numItems;
  int[] lastRow, thisRow;
  BitSet[] pickset;
    
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
        numItems = Integer.parseInt(firstLine[0]);
        capacity = Integer.parseInt(firstLine[1]);
        int[] values = new int[numItems];
        int[] weights = new int[numItems];

        for(int i= 0; i < numItems; i++){
          String line = lines.get(i+1);
          String[] parts = line.split("\\s+");
          values[i] = Integer.parseInt(parts[0]);
          weights[i] = Integer.parseInt(parts[1]);
        }
        
        // dynamic programming ... recycle rows to reduce memory; bitset to track selections
        // note: may need to set -Xmx parameter in java to get more memory/heap

        pickset = new BitSet[numItems + 1];
        for (int i = 0; i <= numItems; i++) { pickset[i] = new BitSet(32); }
        lastRow = new int[capacity + 1];
        thisRow = new int[capacity + 1];
        for (int x = 0; x <= capacity; x++) { lastRow[x] = 0; }
        for (int i = 1; i <= numItems; i++) {
          for (int x = 0; x <= capacity; x++) {
            thisRow[x] = lastRow[x];
            if (x >= weights[i-1]) {
              int z = lastRow[x - weights[i-1]] + values[i-1];
              if (z > thisRow[x]) { 
                thisRow[x] = z; 
                pickset[i].set(x);
              }
            }
          }
          if (i < numItems) for (int x = 0; x <= capacity; x++) { lastRow[x] = thisRow[x]; }
        }
        
        // reconstruct the solution 
        int value = 0;
        int taken[] = new int[numItems];
        int x = capacity;
        for (int i = numItems; i > 0; i--) {
          int nextX = x;
          if (pickset[i].get(x)) {
            taken[i-1] = 1;
            value += values[i-1];
            nextX = x - weights[i-1];            
          }
          x = nextX;
        }
  
      
        // prepare the solution in the specified output format
        System.out.println(value+" 0");
        for(int i=0; i < numItems; i++){
          System.out.print(taken[i]+" ");
        }
        System.out.println(""); 

    }

}