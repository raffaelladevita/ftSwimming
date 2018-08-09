package org.clas;

import java.util.ArrayList;
import javax.swing.JFrame;

import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F; 
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.math.F1D;
import org.jlab.clas.physics.Particle;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.utils.groups.IndexedList;

import cnuphys.rk4.IStopper;
import cnuphys.swim.Swimmer;
import cnuphys.swim.SwimTrajectory;
import cnuphys.magfield.Solenoid;
import java.io.File;
import java.io.FileNotFoundException;


/*
* Analyze tracking results
* TODO:  
* @author devita
*/

public class swimTest {
    
    static final boolean debug=false;
    
    String  analysisName = "Fit Test";
    
    IndexedList<DataGroup> dataGroups      = new IndexedList<DataGroup>(1);
    EmbeddedCanvasTabbed   canvasTabbed    = null;
    EmbeddedCanvas         canvas          = null;
    ArrayList<String>      canvasTabNames  = new ArrayList<String>();
    
    // these correspond to Joseph's two-particle event generater:
    static final int electronSector=1;
    static final int hadronSector=3;
    
    boolean isForwardTagger=false;
    boolean isCentral=false;
    int fdCharge = 0;
    
    int nNegTrackEvents = 0;
    int nTwoTrackEvents = 0;
    
    int nEvents = 0;
    
    double ebeam = 10.6;
    
    static Swimmer solSwimmer=null;
    static Solenoid solMap=null;

    final double swimStepSize=0.0002; // m
    final double swimStepSave=0.0002; // m
    final double swimMaxDistance=2.5; // m
    
    public class FTStopper implements IStopper {
        private double vertexPlane=0; // m
        public boolean stopIntegration(double t,double[] y) {
            final double xx=y[SwimTrajectory.X_IDX];
            final double yy=y[SwimTrajectory.Y_IDX];
            final double zz=y[SwimTrajectory.Z_IDX];
            final double radius = Math.sqrt(xx*xx+yy*yy);
            return zz<=vertexPlane;
        }
        public double getFinalT() { return 0; }
        public void setFinalT(double finalT) { }
    }

    public static void main(String arg[]) throws FileNotFoundException{
        
//  System.setProperty("java.awt.headless", "true"); // this should disable the Xwindow requirement
    GStyle.getAxisAttributesX().setTitleFontSize(24);
    GStyle.getAxisAttributesX().setLabelFontSize(18);  
    GStyle.getAxisAttributesY().setTitleFontSize(24);
    GStyle.getAxisAttributesY().setLabelFontSize(18);
    GStyle.getAxisAttributesZ().setLabelFontSize(14);
//  GStyle.setPalette("kDefault");
    GStyle.getAxisAttributesX().setLabelFontName("Avenir");
    GStyle.getAxisAttributesY().setLabelFontName("Avenir");
    GStyle.getAxisAttributesZ().setLabelFontName("Avenir");
    GStyle.getAxisAttributesX().setTitleFontName("Avenir");
    GStyle.getAxisAttributesY().setTitleFontName("Avenir");
    GStyle.getAxisAttributesZ().setTitleFontName("Avenir");
    GStyle.setGraphicsFrameLineWidth(1);
    GStyle.getH1FAttributes().setLineWidth(1);
    GStyle.getH1FAttributes().setOptStat("1111");
        
    
        swimTest ttest = new swimTest();
       
        File solMapFile=new File("/Users/devita/NetBeansProjects/clas12-offline-software/coatjava/etc/data/magfield/Symm_solenoid_r601_phi1_z1201_13June2018.dat");
        solMap=Solenoid.fromBinaryFile(solMapFile);
        swimTest.solSwimmer=new Swimmer(solMap);
        
        ttest.setAnalysisTabNames("test","FitTest2D");
        ttest.createHistos();
  
        HipoDataSource reader = new HipoDataSource();
        reader.open("/Users/devita/out_out.hipo");

        while (reader.hasEvent()) {
            DataEvent event = reader.getNextEvent();
//            event.show();
        ttest.processEvent(event);
        }
        reader.close();

        JFrame frame = new JFrame("ftSwimming");
        frame.setSize(1200, 800);
        frame.add(ttest.canvasTabbed);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);

        ttest.analyze();
        ttest.plotHistos();
   //     System.out.println("done");    

    }    

    private void processEvent(DataEvent event) {
 
        nEvents++;
        if((nEvents%1000) == 0) System.out.println("Analyzed " + nEvents + " events");
    
        DataBank recRun = null;
        DataBank mcBank = null;
        DataBank ftCalClusters = null;
        if(event.hasBank("RUN::config"))            recRun = event.getBank("RUN::config");
        if(event.hasBank("MC::Particle"))          mcBank  = event.getBank("MC::Particle");
        if(event.hasBank("FTCAL::clusters")) ftCalClusters = event.getBank("FTCAL::clusters");
        int ev = recRun.getInt("event",0);

        if (ftCalClusters != null && mcBank != null) {
            Particle p1 = null;
            Particle p2 = null;
//            mcBank.show();          
                double px = mcBank.getFloat("px", 0); 
                double py = mcBank.getFloat("py", 0); 
                double pz = mcBank.getFloat("pz", 0); 
// Magnitude/Normalization
                double length = Math.sqrt(px*px+py*py+pz*pz);   
                double p_x = px/length; 
                double p_y = py/length;
                double p_z = pz/length; 
  //Angle                      
                double phi = Math.toDegrees(Math.atan2(p_y, p_x));              
                double theta = Math.toDegrees(Math.acos(p_z)); 
// double lengthG = Math.sqrt(px*px+py*py+pz*pz);   
               
        for(int i=0; i<ftCalClusters.rows(); i++) {             
                                
                double energy  = ftCalClusters.getFloat("energy", i);
                double energyR = ftCalClusters.getFloat("recEnergy", i);
                int    size    = ftCalClusters.getInt("size", i);
                double x       = ftCalClusters.getFloat("x", i);
                double y       = ftCalClusters.getFloat("y", i);
                double z       = ftCalClusters.getFloat("z", i);
                double path = Math.sqrt(x*x+y*y+z*z);
                double time = ftCalClusters.getFloat("time", i)-path/29.97;
                double cx=x/path;
                double cy=y/path;
                double cz=z/path;

                double datagroupphi = Math.toDegrees(Math.atan2(cy, cx)); 
                double datagrouptheta = Math.toDegrees(Math.acos(cz)); 
                
                double swimphi   = Math.toDegrees(Math.atan2(-cy, -cx)); 
                double swimtheta = Math.toDegrees(Math.acos(-cz)); 
                
                double pathG = Math.sqrt(x*x+y*y+z*z);
               
                
                
                FTStopper ftStopper = new FTStopper();
                
                int charge=-1;

                SwimTrajectory traj=solSwimmer.swim(charge,x/100,y/100,z/100,
                        energy,swimtheta,swimphi,
                        ftStopper,swimMaxDistance,swimStepSize,swimStepSave);
                traj.computeBDL(solMap);
//                for(int istep=0; istep<traj.size(); istep++) {
//                    System.out.println(traj.get(istep)[SwimTrajectory.Z_IDX] + " " + traj.get(istep)[SwimTrajectory.PATHLEN_IDX] );
//                }
                double[] lastStep = traj.get(traj.size()-1);
                double pathLength = lastStep[SwimTrajectory.PATHLEN_IDX]*100; //cm
                double zvertex    = lastStep[SwimTrajectory.Z_IDX]*100;
                double vertexCX = -lastStep[SwimTrajectory.DIRCOSX_IDX];
                double vertexCY = -lastStep[SwimTrajectory.DIRCOSY_IDX];
                double vertexCZ = -lastStep[SwimTrajectory.DIRCOSZ_IDX];
                double vertexPhi   = Math.toDegrees(Math.atan2(vertexCY, vertexCX)); 
                double vertexTheta = Math.toDegrees(Math.acos(vertexCZ)); 

                double pathLengthG = lastStep[SwimTrajectory.PATHLEN_IDX]*100; //cm
              
                
//                System.out.println(phi + " " + " " + datagroupphi +" "+ vertexPhi);

                                               
                             
                dataGroups.getItem(1).getH2F("hi_cal_E_phi").fill(energy,datagroupphi);
                dataGroups.getItem(1).getH2F("hi_cal_E_theta").fill(energy,datagrouptheta); 
                dataGroups.getItem(1).getH2F("Diff #phi Vs Energy").fill(energy, -datagroupphi +phi); 
                dataGroups.getItem(1).getH2F("Swim: Diff #phi Vs Energy").fill(energy, -vertexPhi + phi); 
                dataGroups.getItem(1).getH2F("Diff #theta Vs Energy").fill(energy,-datagrouptheta +theta);
//                        System.out.println(energy + " " + Math.toDegrees(Math.acos(cz) - Math.acos(p_z)));             
                dataGroups.getItem(1).getH1F("hi_cal_e_ch").fill(energy);
                dataGroups.getItem(1).getH1F("hi_cal_theta_ch").fill(Math.toDegrees(Math.acos(cz)));
                dataGroups.getItem(1).getH1F("hi_cal_phi_ch").fill(datagroupphi);                 
dataGroups.getItem(1).getH2F ("Vertix Pathlength vs Energy").fill(energy, +pathLengthG -pathG);
dataGroups.getItem(1).getH2F ("Swim: mctheta-vtheta Vs Energy").fill(energy, vertexTheta - theta);

                
        }
    }
}
    private void createHistos() {
        // Calorimeter
        DataGroup dc_calo = new DataGroup(2,20);
       
        H2F histogram2d = new H2F("hi_cal_E_theta", 100, 0, 5, 100, 2, 5); 
        H2F histogram2d2 = new H2F("hi_cal_E_theta2", 100, 0, 7, 100, 2, 5.5);  
        H2F histogramh2d = new H2F("hi_cal_E_phi", 100, 0, 5, 100, 0, 190);           
        H2F histogramh2d2 = new H2F("Diff #phi Vs Energy", 100, 1.0, 4.9, 100,3.0, 25);  
        H2F histogramd2D = new H2F("Diff #theta Vs Energy", 100, 0.5, 4.9, 100,-0.4, 0.5);  
        H2F histogramhSwim2d2 = new H2F("Swim: Diff #phi Vs Energy", 100, 1.0, 4.9, 100,-5, 5);  
        H2F histogramPath = new H2F("Vertix Pathlength vs Energy", 100, 0.5, 5, 100, -0.03, 0.10); 
        H2F histogramhSwim2D = new H2F("Swim: mctheta-vtheta Vs Energy", 100, 0.5, 4.9, 100,-0.4, 0.5); 

        
        histogram2d.setTitle("hi_cal_E_theta2D"); 
        histogram2d.setTitleX("Energy"); 
        histogram2d.setTitleY("Angle");	
        histogram2d2.setTitle("hi_cal_E_theta2"); 
        histogram2d2.setTitleX("Energy");
        histogram2d2.setTitleY("Angle");
        histogramh2d.setTitle("#phi 2D");
        histogramh2d.setTitleX("#phi");
        histogramh2d.setTitleY("Angle"); 
        histogramh2d2.setTitle("Mc #phi -rec #phi Vs. Energy");           
        histogramh2d2.setTitleX("Energy");
        histogramh2d2.setTitleY("Difference #phi");  
        histogramd2D.setTitle("Mc #theta -rec #theta Vs. Energy");        
        histogramd2D.setTitleX("Energy");
        histogramd2D.setTitleY("Difference #theta"); 
        histogramhSwim2d2.setTitle("Swim: Diff #phi Vs Energy");     
        histogramhSwim2d2.setTitleX("Energy"); 
        histogramhSwim2d2.setTitleY("Diff #phi"); 
        histogramhSwim2D.setTitle("Swim: Diff #theta Vs Energy");     
        histogramhSwim2D.setTitleX("Energy"); 
        histogramhSwim2D.setTitleY("mctheta-vtheta"); 
        histogramPath.setTitle("Vertix Pathlength vs Energy");     
        histogramPath.setTitleX("Energy"); 
        histogramPath.setTitleY("vertixpathlength"); 

                
                double phicorr0 = 4.918;
                double phicorr1 = -3.828; 
                double phicorr2 = 3.841; 
                double phicorr3 = -1.256; 
                double phicorr4 = 2.874; 
                double phicorr5 = -0.2195; 
                double thetacorr0 = 1.797; 
                double thetacorr1 = -4.485;
                double thetacorr2 = -0.8671; 
                double thetacorr3 = -1.078;  
        
    F1D func_theta = new F1D("func_theta","exp([t0] + [t1] * x)+exp([t2] + [t3] * x)", 1.0, 4.7); 
        func_theta.setParameter(0, thetacorr0);
        func_theta.setParameter(1, thetacorr1);
        func_theta.setParameter(2, thetacorr2);
        func_theta.setParameter(3, thetacorr3);
        func_theta.setLineWidth(2);
        
//   F1D func_thetaA = new F1D("func_theta","exp([t0] + [t1] * x)+exp([t2] + [t3] * x)", 1.0, 4.7); 
//        func_thetaA.setParameter(0, thetacorr0);
//        func_thetaA.setParameter(1, thetacorr1);
//        func_thetaA.setParameter(2, thetacorr2);
//        func_thetaA.setParameter(3, thetacorr3);
//        func_thetaA.setLineWidth(2);
//        func_thetaA.setLineColor(4);     
//        ftheta.setFillColor(2);
  
    F1D func_phi = new F1D("func_phi", "exp([p0] + [p1] * x)+exp([p2] + [p3] * x)+exp([p4] + [p5] * x)", 1.3, 4.7);
        func_phi.setParameter(0, phicorr0);
        func_phi.setParameter(1, phicorr1);
        func_phi.setParameter(2, phicorr2);
        func_phi.setParameter(3, phicorr3);
        func_phi.setParameter(4, phicorr4);
        func_phi.setParameter(5, phicorr5);
        func_phi.setLineWidth(2);
        
        H1F hi_cal_e_ch = new H1F("hi_cal_e_ch", "E (GeV/c^2)", "Counts", 100, 0, 5); 
        hi_cal_e_ch.setFillColor(2);
        H1F hi_cal_nclusters = new H1F("hi_cal_nclusters", 5, 0, 5.5);    
        hi_cal_nclusters.setFillColor(44);
        H1F hi_cal_theta_ch = new H1F("hi_cal_theta_ch","#theta (deg)", "Counts", 100, 2, 5); 
        hi_cal_theta_ch.setFillColor(2);
        H1F hi_cal_phi_ch = new H1F("hi_cal_phi_ch", "#phi (deg)", "Counts", 100, -200,200); 
        hi_cal_phi_ch.setFillColor(2);              
        
        dc_calo.addDataSet(func_phi,         0);
        dc_calo.addDataSet(func_theta,       1);
        dc_calo.addDataSet(histogram2d,      2);
        dc_calo.addDataSet(hi_cal_e_ch,      3);
        dc_calo.addDataSet(hi_cal_theta_ch,  4);
        dc_calo.addDataSet(hi_cal_phi_ch,    5);
        dc_calo.addDataSet(histogram2d2,     6);
        dc_calo.addDataSet(histogramh2d,     7); 
        dc_calo.addDataSet(histogramh2d2,    8);         
        dc_calo.addDataSet(histogramd2D,     9);
        dc_calo.addDataSet(histogramhSwim2d2,10);    
        dc_calo.addDataSet(histogramPath,    11);
        dc_calo.addDataSet(histogramhSwim2D, 12);
//      dc_calo.addDataSet(func_thetaA,      13);
     
        dataGroups.add(dc_calo, 1);
    }
    private void plotHistos() {

        canvasTabbed.getCanvas("test").divide(3,2);
        canvasTabbed.getCanvas("test").setGridX(false);
        canvasTabbed.getCanvas("test").setGridY(false); 
        canvasTabbed.getCanvas("FitTest2D").divide(2,2);
        canvasTabbed.getCanvas("FitTest2D").setGridX(false);
        canvasTabbed.getCanvas("FitTest2D").setGridY(false);
        canvasTabbed.getCanvas("test").cd(0);
        canvasTabbed.getCanvas("test").draw(dataGroups.getItem(1).getH1F("hi_cal_phi_ch"));
        canvasTabbed.getCanvas("test").cd(1);
        canvasTabbed.getCanvas("test").draw(dataGroups.getItem(1).getH2F("hi_cal_E_phi"));
        canvasTabbed.getCanvas("test").cd(2);
        canvasTabbed.getCanvas("test").draw(dataGroups.getItem(1).getH2F("Diff #phi Vs Energy"));
        canvasTabbed.getCanvas("test").draw(dataGroups.getItem(1).getF1D("func_phi"),"same");  
        canvasTabbed.getCanvas("test").cd(3);
        canvasTabbed.getCanvas("test").draw(dataGroups.getItem(1).getH2F("Vertix Pathlength vs Energy"));
        canvasTabbed.getCanvas("test").cd(4);
        canvasTabbed.getCanvas("test").draw(dataGroups.getItem(1).getH2F("Swim: mctheta-vtheta Vs Energy"));
        canvasTabbed.getCanvas("test").draw(dataGroups.getItem(1).getF1D("func_theta"),"same");  
        canvasTabbed.getCanvas("test").cd(5); 
        canvasTabbed.getCanvas("test").draw(dataGroups.getItem(1).getH2F("Swim: Diff #phi Vs Energy"));

        canvasTabbed.getCanvas("FitTest2D").cd(0);
        canvasTabbed.getCanvas("FitTest2D").draw(dataGroups.getItem(1).getH2F("hi_cal_E_theta"));
        canvasTabbed.getCanvas("FitTest2D").cd(1);
        canvasTabbed.getCanvas("FitTest2D").draw(dataGroups.getItem(1).getH1F("hi_cal_e_ch"));
        canvasTabbed.getCanvas("FitTest2D").cd(2);
        canvasTabbed.getCanvas("FitTest2D").draw(dataGroups.getItem(1).getH1F("hi_cal_theta_ch"));
        canvasTabbed.getCanvas("FitTest2D").cd(3);
//canvasTabbed.getCanvas("FitTest2D").draw(dataGroups.getItem(1).getF1D("func_thetaA"),"same");  
//canvasTabbed.getCanvas("FitTest2D").draw(dataGroups.getItem(1).getH2F("Swim: mctheta-vtheta Vs Energy"));
        canvasTabbed.getCanvas("FitTest2D").draw(dataGroups.getItem(1).getH2F("Diff #theta Vs Energy"));
        canvasTabbed.getCanvas("FitTest2D").draw(dataGroups.getItem(1).getF1D("func_theta"),"same");
        canvasTabbed.getCanvas("FitTest2D").cd(4); 
    } 
    private void analyze() {
//        System.out.println("Updating TBT");
//        //fitting negative tracks vertex
//        this.fitVertex(dataGroups.getItem(1).getH1F("hi_vz_neg_cut"), dataGroups.getItem(1).getF1D("f1_vz_neg"));
//        //fitting positive tracks vertex
//        this.fitVertex(dataGroups.getItem(2).getH1F("hi_vz_pos_cut"), dataGroups.getItem(2).getF1D("f1_vz_pos"));
//        // fitting MC comparisons
//        this.fitMC(dataGroups.getItem(3).getH1F("hi_dp_pos"),     dataGroups.getItem(3).getF1D("f1_dp_pos"));
//        this.fitMC(dataGroups.getItem(3).getH1F("hi_dtheta_pos"), dataGroups.getItem(3).getF1D("f1_dtheta_pos"));
//        this.fitMC(dataGroups.getItem(3).getH1F("hi_dphi_pos"),   dataGroups.getItem(3).getF1D("f1_dphi_pos"));
//        this.fitMC(dataGroups.getItem(3).getH1F("hi_dvz_pos"),    dataGroups.getItem(3).getF1D("f1_dvz_pos"));     
//        this.fitMC(dataGroups.getItem(3).getH1F("hi_dp_neg"),     dataGroups.getItem(3).getF1D("f1_dp_neg"));
//        this.fitMC(dataGroups.getItem(3).getH1F("hi_dtheta_neg"), dataGroups.getItem(3).getF1D("f1_dtheta_neg"));
//        this.fitMC(dataGroups.getItem(3).getH1F("hi_dphi_neg"),   dataGroups.getItem(3).getF1D("f1_dphi_neg"));
//        this.fitMC(dataGroups.getItem(3).getH1F("hi_dvz_neg"),    dataGroups.getItem(3).getF1D("f1_dvz_neg"));     
//        this.fitMCSlice(dataGroups.getItem(3).getH2F("hi_dp_p_neg"),dataGroups.getItem(3).getGraph("gr_dp_p_neg"));
//        this.fitMCSlice(dataGroups.getItem(3).getH2F("hi_dp_theta_neg"),dataGroups.getItem(3).getGraph("gr_dp_theta_neg"));
//        this.fitMCSlice(dataGroups.getItem(3).getH2F("hi_dp_phi_neg"),dataGroups.getItem(3).getGraph("gr_dp_phi_neg"));
    }
    private void fitVertex(H1F hivz, F1D f1vz) {
        double mean  = hivz.getDataX(hivz.getMaximumBin());
        double amp   = hivz.getBinContent(hivz.getMaximumBin());
        double sigma = 1.;
        if(hivz.getEntries()>500) { // first fits 
            sigma = Math.abs(f1vz.getParameter(2));       
        }
        f1vz.setParameter(0, amp);
        f1vz.setParameter(1, mean);
        f1vz.setParameter(2, sigma);
        f1vz.setRange(mean-2.*sigma,mean+2.*sigma);
        DataFitter.fit(f1vz, hivz, "FitTest2D"); //No options uses error for sigma 
        hivz.setFunction(null);
    }

    private void fitMC(H1F himc, F1D f1mc) {
        double mean  = himc.getDataX(himc.getMaximumBin());
        double amp   = himc.getBinContent(himc.getMaximumBin());
        double sigma = himc.getRMS()/2;
        f1mc.setParameter(0, amp);
        f1mc.setParameter(1, mean);
        f1mc.setParameter(2, sigma);
        f1mc.setRange(mean-2.*sigma,mean+2.*sigma);
        DataFitter.fit(f1mc, himc, "FitTest2D"); //No options uses error for sigma 
        sigma = Math.abs(f1mc.getParameter(2));  
        f1mc.setRange(mean-2.*sigma,mean+2.*sigma);
        DataFitter.fit(f1mc, himc, "FitTest2D"); //No options uses error for sigma 
        himc.setFunction(null);
    }
    
    private void fitMCSlice(H2F himc, GraphErrors grmc) {
        grmc.reset();
        ArrayList<H1F> hslice = himc.getSlicesX();
        for(int i=0; i<hslice.size(); i++) {
            double  x = himc.getXAxis().getBinCenter(i);
            double ex = 0;
            double  y = hslice.get(i).getRMS();
            double ey = 0;
            double mean  = hslice.get(i).getDataX(hslice.get(i).getMaximumBin());
            double amp   = hslice.get(i).getBinContent(hslice.get(i).getMaximumBin());
            double sigma = hslice.get(i).getRMS()/2;
            F1D f1 = new F1D("f1slice","[amp]*gaus(x,[mean],[sigma])", hslice.get(i).getDataX(0), hslice.get(i).getDataX(hslice.get(i).getDataSize(1)-1));
            f1.setParameter(0, amp);
            f1.setParameter(1, mean);
            f1.setParameter(2, sigma);
            f1.setRange(mean-2.*sigma,mean+2.*sigma);
            DataFitter.fit(f1, hslice.get(i), "FitTest2D"); //No options uses error for sigma 
            if(amp>50) grmc.addPoint(x, f1.getParameter(2), ex, f1.parameter(2).error());
        }

    }  
    public void setAnalysisTabNames(String... names) {
        for(String name : names) {
            canvasTabNames.add(name);
        }
            canvasTabbed = new EmbeddedCanvasTabbed(names);
    }    
    public void printCanvas(String dir, String name) {
        // print canvas to files
        for(int tab=0; tab<canvasTabNames.size(); tab++) {
            String fileName = dir + "/" + this.analysisName + "_" + name + "." + tab + ".png";
            System.out.println(fileName);
            canvasTabbed.getCanvas(canvasTabNames.get(tab)).save(fileName);
        }
    }
    public boolean testMCpart(Particle mc, Particle rec) {
        if(Math.abs(mc.px()-rec.px())<2.5 &&
           Math.abs(mc.py()-rec.py())<2.5 &&
           Math.abs(mc.pz()-rec.pz())<2.5) return true;
        else return false;
    }
}