/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas;

import java.util.ArrayList;
import javax.swing.JFrame;
import org.jlab.clas.physics.LorentzVector;

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
import org.jlab.groot.group.DataGroup;
import org.jlab.utils.groups.IndexedList;

/**
 *
 * Analyze tracking results
 *
 * TODO:  
 * 
 * @author devita
 */
/**
 *
 * @author devita
 */

public class swimTest {

    static final boolean debug=false;
    
    String  analysisName = "FT";
    
    IndexedList<DataGroup> dataGroups      = new IndexedList<DataGroup>(1);
    EmbeddedCanvasTabbed   canvasTabbed    = null;
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


    public static void main(String arg[]){
        
//        System.setProperty("java.awt.headless", "true"); // this should disable the Xwindow requirement
        GStyle.getAxisAttributesX().setTitleFontSize(24);
        GStyle.getAxisAttributesX().setLabelFontSize(18);
        GStyle.getAxisAttributesY().setTitleFontSize(24);
        GStyle.getAxisAttributesY().setLabelFontSize(18);
        GStyle.getAxisAttributesZ().setLabelFontSize(14);
//        GStyle.setPalette("kDefault");
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
        
        ttest.setAnalysisTabNames("FT");
        ttest.createHistos();
        


        HipoDataSource reader = new HipoDataSource();
        reader.open("/Users/devita/out_out.hipo");

        while (reader.hasEvent()) {
            DataEvent event = reader.getNextEvent();
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
        

    }


    private void processEvent(DataEvent event) {

        nEvents++;
        if((nEvents%10000) == 0) System.out.println("Analyzed " + nEvents + " events");
    
        DataBank recRun    = null;
        DataBank mcBank    = null;
        DataBank ftCalClusters = null;
        if(event.hasBank("RUN::config"))            recRun      = event.getBank("RUN::config");
        if(event.hasBank("MC::Particle"))          mcBank       = event.getBank("MC::Particle");
        if(event.hasBank("FTCAL::clusters"))      ftCalClusters = event.getBank("FTCAL::clusters");
        int ev = recRun.getInt("event",0);

        if (ftCalClusters != null && mcBank != null) {
            Particle p1 = null;
            Particle p2 = null;
            dataGroups.getItem(1).getH1F("hi_cal_nclusters").fill(ftCalClusters.rows());
            for(int i=0; i<ftCalClusters.rows(); i++) {
                double energy  = ftCalClusters.getFloat("energy", i);
                double energyR  = ftCalClusters.getFloat("recEnergy", i);
                int    size     = ftCalClusters.getInt("size", i);
                double x = ftCalClusters.getFloat("x", i);
                double y = ftCalClusters.getFloat("y", i);
                double z = ftCalClusters.getFloat("z", i);
                double path     = Math.sqrt(x*x+y*y+z*z);
                double time     = ftCalClusters.getFloat("time", i)-path/29.97;
                double cx=x/path;
                double cy=y/path;
                double cz=z/path;
                dataGroups.getItem(1).getH1F("hi_cal_e_ch").fill(energy);
                dataGroups.getItem(1).getH1F("hi_cal_theta_ch").fill(Math.toDegrees(Math.acos(cz)));
                dataGroups.getItem(1).getH1F("hi_cal_phi_ch").fill(Math.toDegrees(Math.atan2(cy,cx)));
               
            }
        }
    }

    
    private void createHistos() {
    
        // Calorimeter
        DataGroup dc_calo = new DataGroup(2,2);
        H1F hi_cal_nclusters = new H1F("hi_cal_nclusters", "N. Clusters", "Counts", 5, 0, 5);    
        hi_cal_nclusters.setFillColor(44);
        H1F hi_cal_e_ch = new H1F("hi_cal_e_ch", "E (GeV)", "Counts", 100, 0, 12); 
        hi_cal_e_ch.setFillColor(2);
        H1F hi_cal_theta_ch = new H1F("hi_cal_theta_ch","#theta (deg)", "Counts", 100, 2,  6); 
        hi_cal_theta_ch.setFillColor(2);
        H1F hi_cal_phi_ch = new H1F("hi_cal_phi_ch", "#phi (deg)", "Counts", 100, -180,180); 
        hi_cal_phi_ch.setFillColor(2);
        dc_calo.addDataSet(hi_cal_nclusters, 0);
        dc_calo.addDataSet(hi_cal_e_ch,      1);
        dc_calo.addDataSet(hi_cal_theta_ch,  2);
        dc_calo.addDataSet(hi_cal_phi_ch,    3);
        dataGroups.add(dc_calo, 1);


    }

    private void plotHistos() {

        canvasTabbed.getCanvas("FT").divide(2,2);
        canvasTabbed.getCanvas("FT").setGridX(false);
        canvasTabbed.getCanvas("FT").setGridY(false);
        canvasTabbed.getCanvas("FT").cd(0);
        canvasTabbed.getCanvas("FT").draw(dataGroups.getItem(1).getH1F("hi_cal_nclusters"));
        canvasTabbed.getCanvas("FT").cd(1);
        canvasTabbed.getCanvas("FT").draw(dataGroups.getItem(1).getH1F("hi_cal_e_ch"));
        canvasTabbed.getCanvas("FT").cd(2);
        canvasTabbed.getCanvas("FT").draw(dataGroups.getItem(1).getH1F("hi_cal_theta_ch"));
        canvasTabbed.getCanvas("FT").cd(3);
        canvasTabbed.getCanvas("FT").draw(dataGroups.getItem(1).getH1F("hi_cal_phi_ch"));

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
        DataFitter.fit(f1vz, hivz, "Q"); //No options uses error for sigma 
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
        DataFitter.fit(f1mc, himc, "Q"); //No options uses error for sigma 
        sigma = Math.abs(f1mc.getParameter(2));  
        f1mc.setRange(mean-2.*sigma,mean+2.*sigma);
        DataFitter.fit(f1mc, himc, "Q"); //No options uses error for sigma 
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
            DataFitter.fit(f1, hslice.get(i), "Q"); //No options uses error for sigma 
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
