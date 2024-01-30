package ru.biosoft.physicell.core;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.ui.Visualizer;
import ru.biosoft.physicell.ui.Visualizer.Section;

public class Model
{
    private List<Visualizer> visualizers = new ArrayList<Visualizer>();
    private Microenvironment m;
    private Map<String, String> parameters = new HashMap<>();
    private double tMax;
    private double diffusion_dt;
    private double mechanics_dt = 0.1;
    private double phenotype_dt = 6.0;
    private double full_save_interval;
    private String resultFolder;
    private boolean enableFullSaves;

    public Iterable<Visualizer> getVisualizers()
    {
        return visualizers;
    }

    public void setResultFolder(String folder)
    {
        this.resultFolder = folder;
    }

    public Visualizer addVisualizer(int zSlice, String name)
    {
        Visualizer visualizer = new Visualizer( resultFolder, name, Section.Z, zSlice );
        //        visualizer.setDrawDensity( false );
        visualizer.setSaveImage( false );
        //        visualizer.setColorPhase( "Ki67-", Color.lightGray );
        //        visualizer.setColorPhase( "Ki67+ (premitotic)", Color.green );
        //        visualizer.setColorPhase( "Ki67+ (postmitotic)", new Color( 0, 128, 0 ) );
        //        visualizer.setColorPhase( "Apoptotic", Color.red );
        this.visualizers.add( visualizer );
        return visualizer;
    }

    public Model()
    {
        m = new Microenvironment();
    }

    public Model(Microenvironment m)
    {
        this.m = m;
    }

    public Microenvironment getMicroenvironment()
    {
        return m;
    }

    public void createContainer(double voxelSize)
    {
        CellContainer.createCellContainer( m, voxelSize );
    }

    public void simulate() throws Exception
    {
        double curTime = 0;
        double next_full_save_time = 0;
        boolean enable_legacy_saves = false;
        //        boolean enable_full_saves = false;
        File reportFile;
        int full_output_index = 0;

        for( Visualizer listener : visualizers )
            listener.init();

        try
        {
            while( curTime < tMax + 0.1 * diffusion_dt )
            {
                // save data if it's time. 
                if( Math.abs( curTime - next_full_save_time ) < 0.01 * diffusion_dt )
                {
                    //                    display_simulation_status( std::cout ); 
                    if( enable_legacy_saves )
                    {
                        //                                    log_output( curTime, full_output_index, m, reportFile);
                    }

                    if( enableFullSaves )
                    {
                        for( Visualizer listener : visualizers )
                            listener.saveResult( m, curTime );
                        //                                    sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 

                        //                                    save_PhysiCell_to_MultiCellDS_v2( filename , m , curTime ); 
                    }

                    full_output_index++;
                    next_full_save_time += full_save_interval;
                }
                //                
                //                // save SVG plot if it's time
                //                if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
                //                {
                //                    if( PhysiCell_settings.enable_SVG_saves == true )
                //                    {   
                //                        sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
                //                        SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
                //                        
                //                        PhysiCell_globals.SVG_output_index++; 
                //                        PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
                //                    }
                //                }
                //
                // update the microenvironment
                m.simulate_diffusion_decay( diffusion_dt );

                // run PhysiCell 
                //                double diffusion_dt = 0.01;
                //                double mechanics_dt = 0.1;
                //                double phenotype_dt = 6.0;
                ( (CellContainer)m.agentContainer ).updateAllCells( m, curTime, phenotype_dt, mechanics_dt, diffusion_dt );

                /*
                  Custom add-ons could potentially go here. 
                */
                curTime += diffusion_dt;
            }

            for( Visualizer listener : visualizers )
                listener.finish();

            if( enable_legacy_saves )
            {
                //                            log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
                //                            report_file.close();
            }
        }

        //        
        //        // save a final simulation snapshot 
        //        
        //        sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
        //        save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
        //        
        //        sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
        //        SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
        //        
        //        // timer 
        //        
        //        std::cout << std::endl << "Total simulation runtime: " << std::endl; 
        //        BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 
        catch( Exception ex )
        {
            ex.printStackTrace();
        }

    }

    public void addParameter(String name, String val)
    {
        this.parameters.put( name, val );
    }

    public String getParameter(String name)
    {
        return parameters.get( name );
    }

    public int getParameterInt(String name)
    {
        return Integer.parseInt( parameters.get( name ) );
    }

    public double getParameterDouble(String name)
    {
        return Double.parseDouble( parameters.get( name ) );
    }

    public void setTMax(double tMax)
    {
        this.tMax = tMax;
    }

    public void setDiffusionDt(double diffusion_dt)
    {
        this.diffusion_dt = diffusion_dt;
    }

    public void setMechanicsDt(double mechanics_dt)
    {
        this.mechanics_dt = mechanics_dt;
    }

    public void setPhenotypeDt(double phenotype_dt)
    {
        this.phenotype_dt = phenotype_dt;
    }

    public void setSaveInterval(double full_save_interval)
    {
        this.full_save_interval = full_save_interval;
    }

    public void setEnableFullSaves(boolean enable)
    {
        enableFullSaves = enable;
    }

    public boolean isEnableFullSaves()
    {
        return enableFullSaves;
    }
}
