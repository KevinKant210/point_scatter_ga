import java.io.*;
import java.util.*;
import java.text.*;





public class PointScatter extends FitnessFunction{
    public PointScatter(){
        name = "PointScatter Problem";
    }

    private double calc_distance(double r1,  double r2, double theta_1, double theta_2){
        double dist = 0;
        
        double cos_sub = theta_2-theta_1;

        double cos_comp = Math.cos(Math.toRadians(cos_sub));

        double compa = Math.pow(r1,2);

        double compb = Math.pow(r2,2);

        double compc = 2 * r1 *r2 * cos_comp;
        
        //finish this function at some point
        dist = compa + compb - 2*r1*r2*cos_comp;

        dist = Math.sqrt(dist);

        

        return dist/100;
    }

    public void doRawFitness(Chromo X){
        //divide up the chromosome into n partitions one for each point
        int n = Parameters.numGenes;

        ArrayList<Double> r_set = new ArrayList<>();
        ArrayList<Double> theta_set = new ArrayList<>();

        String sub;
        
        int half = X.chromo.length()/2;

        
        
        for(int i = 0 ;i < n; i++){
           
             sub = X.chromo.substring(i*9, i*9+9);

             double val = X.getSubGeneValue(sub);
             if(val > 100){
                val = val %100;
             }

             

             r_set.add(val);
            
        }

        
        //should probably split the entire chromosome in half and do one for each half this halfs the amount of time im using since I can do one loop

       
        
        for(int i = 0 ;i < n; i++){
           
             sub = X.chromo.substring(half+(i*9), half + (i*9 + 9));
             double val = X.getSubGeneValue(sub);

             

             if(val > 360){
                val = val%360;
             }

             theta_set.add(val);

            
        }



        
        ArrayList<Double> dist = new ArrayList<>();

        for(int i = 0 ; i < n; i++){
            for(int j =0 ; j < n; j++){
                if(i != j){
                    dist.add(calc_distance(r_set.get(i), r_set.get(j), theta_set.get(i), theta_set.get(j)));
                }
            }
        }

        X.rawFitness = Collections.min(dist);
        
        //need to make choices based on handling of invalid children, should probably be done at the crossover level doing at thius level for simplicity
    }

    public void doPrintGenes(Chromo X, FileWriter output) throws java.io.IOException{

		output.write("Chromo:" +  X.chromo + "\n");
        output.write("Raw Fitness :" + String.valueOf(X.rawFitness) + "\n");
		return;
	}
}
