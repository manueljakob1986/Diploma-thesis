package runner;

import com.relteq.sirius.simulator.Node;
import com.relteq.sirius.simulator.Scenario;
import com.relteq.sirius.simulator.Link;
import com.relteq.sirius.simulator.ObjectFactory;
import com.relteq.sirius.simulator.SiriusErrorLog;
import com.relteq.sirius.simulator.SiriusException;
import java.util.ArrayList;
import java.util.List;

public class Runner {

    public static double[][][] getDensities(double dt, double T, Double[][] sourcedemand, Double[][][] splitRatioMatrix) {
	/*
	* sourceDemand[i][j]: i = time index, j = path index
	* splitRatioMatrix[i][j][k] : i = input link, j = output link, k = path
	* output: densities[t][i][k] : t = time, i = links, k = path
	*
	* */
    	
        String configfilename = null;
        int numTimeSteps = (int)(T / dt);
        int nPaths = splitRatioMatrix[0][0].length;
        double[][][] density = new double[numTimeSteps+1][][];

        configfilename = "C:\\Users\\Logge\\workspace\\MJRunner\\data\\mySplitSymmetric_complex_sirius.xml";
        Thread.currentThread().setContextClassLoader(Runner.class.getClassLoader());
        Scenario S = ObjectFactory.createAndLoadScenario(configfilename);

        if(S==null || SiriusErrorLog.haserror()){
            SiriusErrorLog.printErrorMessage();
        }

        List<com.relteq.sirius.jaxb.Node> jaxbNodes = S.getNetworkList().getNetwork().get(0).getNodeList().getNode();

        int nNodes = jaxbNodes.size();
        List<Node> nodes = new ArrayList<Node>(nNodes);

        for (com.relteq.sirius.jaxb.Node node : jaxbNodes) {
            nodes.add(S.getNodeWithCompositeId(null, node.getId()));
        }

        for (Node node : nodes) {
            Link[] inLinks = node.getInput_link();
            Link[] outLinks = node.getOutput_link();
            int inSize = inLinks.length;
            int outSize = outLinks.length;
            if (inSize < 1 || outSize < 1 || nPaths < 1)
                continue;
            double[][][] splitsMatrix = new double[inSize][outSize][nPaths];
            for (int inIndex = 0; inIndex < inSize; inIndex++) {
                Link inLink = inLinks[inIndex];
                int inId = Integer.parseInt(inLink.getId()) - 1;
                for (int outIndex = 0; outIndex < outSize; outIndex++) {
                    Link outLink = outLinks[outIndex];
                    int outId = Integer.parseInt(outLink.getId()) - 1;
                    for (int pathId = 0; pathId < nPaths; pathId++) {
                        splitsMatrix[inIndex][outIndex][pathId] = splitRatioMatrix[inId][outId][pathId];
                    }
                }
            }
            try {
                node.setSplitRatioMatrix(splitsMatrix);
            } catch (SiriusException se) {
                System.out.println(se);
            }
        }

        try {

            // reference to the on ramp link
            Link theonramp = S.getLinkWithCompositeId(null,"1");

            // initialize the run
            S.initialize_run();

            int iTime = 0;
            while (true) {
            	       	
                double[][] temp = S.getDensityForNetwork(null);
                density[iTime] = new double[temp.length][];
                for(int j=0; j<temp.length; j++) {
                    density[iTime][j] = new double[temp[j].length];
                    for(int k=0; k <temp[j].length; k++)
                        density[iTime][j][k] = temp[j][k];
                }

            	if(iTime>=numTimeSteps)
            		break;
            	
                theonramp.setSourcedemandFromVeh(sourcedemand[iTime]);
                S.advanceNSeconds(dt);
                
                iTime++;
            }

        } catch (SiriusException e1) {
            e1.printStackTrace();
        }
        return density;
    }

    public static void main(String[] args){
    		
//    	double dt=3600.0;
//    	double T=18000.0;
//    	Double[][] sourcedemand = new Double[20][1];
//    	sourcedemand[0][0]=0.5;
//    	sourcedemand[1][0]=0.4;
//    	sourcedemand[2][0]=0.3;
//    	sourcedemand[3][0]=0.2;
//    	sourcedemand[4][0]=0.1;
//    	sourcedemand[5][0]=0.05;
//    	sourcedemand[6][0]=0.0;
//    	sourcedemand[7][0]=0.0;
//    	sourcedemand[8][0]=0.0;
//    	sourcedemand[9][0]=0.0;
//    	sourcedemand[10][0]=0.0;
//    	sourcedemand[11][0]=0.0;
//    	sourcedemand[12][0]=0.0;
//    	sourcedemand[13][0]=0.0;
//    	sourcedemand[14][0]=0.0;
//    	sourcedemand[15][0]=0.0;
//    	sourcedemand[16][0]=0.0;
//    	sourcedemand[17][0]=0.0;
//    	sourcedemand[18][0]=0.0;
//    	sourcedemand[19][0]=0.0;
//    	
//    	Double[][][] splitRatioMatrix = new Double[8][8][1];
//    	splitRatioMatrix[0][0][0]=0.0;
//    	splitRatioMatrix[0][1][0]=1.0;
//    	splitRatioMatrix[0][2][0]=0.0;
//    	splitRatioMatrix[0][3][0]=0.0;
//    	splitRatioMatrix[0][4][0]=0.0;
//    	splitRatioMatrix[0][5][0]=0.0;
//    	splitRatioMatrix[0][6][0]=0.0;
//    	splitRatioMatrix[0][7][0]=0.0;
//    	splitRatioMatrix[1][0][0]=0.0;
//    	splitRatioMatrix[1][1][0]=0.0;
//    	splitRatioMatrix[1][2][0]=1.0;
//    	splitRatioMatrix[1][3][0]=0.0;
//    	splitRatioMatrix[1][4][0]=0.0;
//    	splitRatioMatrix[1][5][0]=0.0;
//    	splitRatioMatrix[1][6][0]=0.0;
//    	splitRatioMatrix[1][7][0]=0.0;
//    	splitRatioMatrix[2][0][0]=0.0;
//    	splitRatioMatrix[2][1][0]=0.0;
//    	splitRatioMatrix[2][2][0]=0.0;
//    	splitRatioMatrix[2][3][0]=1.0;
//    	splitRatioMatrix[2][4][0]=0.0;
//    	splitRatioMatrix[2][5][0]=0.0;
//    	splitRatioMatrix[2][6][0]=0.0;
//    	splitRatioMatrix[2][7][0]=0.0;
//    	splitRatioMatrix[3][0][0]=0.0;
//    	splitRatioMatrix[3][1][0]=0.0;
//    	splitRatioMatrix[3][2][0]=0.0;
//    	splitRatioMatrix[3][3][0]=0.0;
//    	splitRatioMatrix[3][4][0]=0.0;
//    	splitRatioMatrix[3][5][0]=0.0;
//    	splitRatioMatrix[3][6][0]=0.0;
//    	splitRatioMatrix[3][7][0]=1.0;
//    	splitRatioMatrix[4][0][0]=0.0;
//    	splitRatioMatrix[4][1][0]=0.0;
//    	splitRatioMatrix[4][2][0]=0.0;
//    	splitRatioMatrix[4][3][0]=0.0;
//    	splitRatioMatrix[4][4][0]=0.0;
//    	splitRatioMatrix[4][5][0]=0.0;
//    	splitRatioMatrix[4][6][0]=0.0;
//    	splitRatioMatrix[4][7][0]=0.0;
//    	splitRatioMatrix[5][0][0]=0.0;
//    	splitRatioMatrix[5][1][0]=0.0;
//    	splitRatioMatrix[5][2][0]=0.0;
//    	splitRatioMatrix[5][3][0]=0.0;
//    	splitRatioMatrix[5][4][0]=0.0;
//    	splitRatioMatrix[5][5][0]=0.0;
//    	splitRatioMatrix[5][6][0]=0.0;
//    	splitRatioMatrix[5][7][0]=0.0;
//    	splitRatioMatrix[6][0][0]=0.0;
//    	splitRatioMatrix[6][1][0]=0.0;
//    	splitRatioMatrix[6][2][0]=0.0;
//    	splitRatioMatrix[6][3][0]=0.0;
//    	splitRatioMatrix[6][4][0]=0.0;
//    	splitRatioMatrix[6][5][0]=0.0;
//    	splitRatioMatrix[6][6][0]=0.0;
//    	splitRatioMatrix[6][7][0]=0.0;
//    	splitRatioMatrix[7][0][0]=0.0;
//    	splitRatioMatrix[7][1][0]=0.0;
//    	splitRatioMatrix[7][2][0]=0.0;
//    	splitRatioMatrix[7][3][0]=0.0;
//    	splitRatioMatrix[7][4][0]=0.0;
//    	splitRatioMatrix[7][5][0]=0.0;
//    	splitRatioMatrix[7][6][0]=0.0;
//    	splitRatioMatrix[7][7][0]=0.0;
//    	
//    	double[][][] x_veh = getDensities(dt,T,sourcedemand,splitRatioMatrix);
//    	
//    	for(int i=0;i<x_veh.length;i++)
//        	for(int j=0;j<x_veh[i].length;j++)
//            	for(int k=0;k<x_veh[i][j].length;k++)
//            		System.out.println("x_veh["+i+"]["+j+"]["+k+"]="+x_veh[i][j][k]);
	
    }

}