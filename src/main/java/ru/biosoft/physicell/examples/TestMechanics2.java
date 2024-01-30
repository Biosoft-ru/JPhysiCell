package ru.biosoft.physicell.examples;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import ru.biosoft.physicell.biofvm.BasicAgent;
import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellContainer;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.StandardModels;
import ru.biosoft.physicell.ui.Visualizer;
import ru.biosoft.physicell.ui.Visualizer.Section;

public class TestMechanics2
{
    /*
    ###############################################################################
    # If you use PhysiCell in your project, please cite PhysiCell and the version #
    # number, such as below:                                                      #
    #                                                                             #
    # We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
    #                                                                             #
    # [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
    #     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
    #     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
    #     DOI: 10.1371/journal.pcbi.1005991                                       #
    #                                                                             #
    # See VERSION.txt or call get_PhysiCell_version() to get the current version  #
    #     x.y.z. Call display_citations() to get detailed information on all cite-#
    #     able software used in your PhysiCell application.                       #
    #                                                                             #
    # Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
    #     as below:                                                               #
    #                                                                             #
    # We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
    # with BioFVM [2] to solve the transport equations.                           #
    #                                                                             #
    # [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
    #     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
    #     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
    #     DOI: 10.1371/journal.pcbi.1005991                                       #
    #                                                                             #
    # [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
    #     llelized diffusive transport solver for 3-D biological simulations,     #
    #     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
    #                                                                             #
    ###############################################################################
    #                                                                             #
    # BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
    #                                                                             #
    # Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
    # All rights reserved.                                                        #
    #                                                                             #
    # Redistribution and use in source and binary forms, with or without          #
    # modification, are permitted provided that the following conditions are met: #
    #                                                                             #
    # 1. Redistributions of source code must retain the above copyright notice,   #
    # this list of conditions and the following disclaimer.                       #
    #                                                                             #
    # 2. Redistributions in binary form must reproduce the above copyright        #
    # notice, this list of conditions and the following disclaimer in the         #
    # documentation and/or other materials provided with the distribution.        #
    #                                                                             #
    # 3. Neither the name of the copyright holder nor the names of its            #
    # contributors may be used to endorse or promote products derived from this   #
    # software without specific prior written permission.                         #
    #                                                                             #
    # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
    # AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
    # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
    # ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
    # LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
    # CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
    # SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
    # INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
    # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
    # ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
    # POSSIBILITY OF SUCH DAMAGE.                                                 #
    #                                                                             #
    ###############################################################################
    */
    double min_voxel_size = 30;
    private static DecimalFormat format = new DecimalFormat( "##.##" );
    private static String resultPath = "C:/Users/Damag/BIOFVM/Mechanics2/Even more/";

    private static boolean outputReport = false;

    private static double[] tumorCenter = new double[] {500, 500, 500};
    static double size = 1000;
    static int zSlice = 500;
    static double dt = 1;
    static double outputInterval = 10.0; // 1.0;
    static double tMax = 20000;

    static double sampleСellRadius = 10;
    static double volume = 2.4943e+03;
    static double sphereRadius = 80;

    private static List<Visualizer> resultListeners = new ArrayList<>();

    static
    {
        Visualizer visualizer = new Visualizer( resultPath, "Z500", Section.Z, zSlice );
        visualizer.setDrawDensity( false );
        visualizer.setSaveImage( false );
        visualizer.setColorPhase( "Ki67-", Color.lightGray );
        visualizer.setColorPhase( "Ki67+ (premitotic)", Color.green );
        visualizer.setColorPhase( "Ki67+ (postmitotic)", new Color( 0, 128, 0 ) );
        visualizer.setColorPhase( "Apoptotic", Color.red );
        resultListeners.add( visualizer );

        Visualizer visualizer2 = new Visualizer( resultPath, "Z590", Section.Z, zSlice - 10 );
        visualizer2.setDrawDensity( false );
        visualizer2.setSaveImage( false );
        visualizer2.setColorPhase( "Ki67-", Color.gray );
        visualizer2.setColorPhase( "Ki67+ (premitotic)", Color.green );
        visualizer2.setColorPhase( "Ki67+ (postmitotic)", new Color( 0, 128, 0 ) );
        visualizer2.setColorPhase( "Apoptotic", Color.red );
        resultListeners.add( visualizer2 );

        Visualizer visualizer3 = new Visualizer( resultPath, "Y500", Section.Y, zSlice );
        visualizer3.setDrawDensity( false );
        visualizer3.setSaveImage( false );
        visualizer3.setColorPhase( "Ki67-", Color.gray );
        visualizer3.setColorPhase( "Ki67+ (premitotic)", Color.green );
        visualizer3.setColorPhase( "Ki67+ (postmitotic)", new Color( 0, 128, 0 ) );
        visualizer3.setColorPhase( "Apoptotic", Color.red );
        resultListeners.add( visualizer3 );
    }

    public static void main(String ... argv) throws Exception
    {
        Microenvironment m = new Microenvironment( "substrate scale", size, 20, "minutes", "microns" );
        CellContainer.createCellContainer( m, 30 );
        m.setDensity( 0, "oxygen", "mmHg", 0, 0 );

        CellDefinition cd = StandardModels.createFromDefault( "tumor cell", 0, m );
        CellDefinition.registerCellDefinition( cd );
        cd.phenotype.cycle = StandardModels.Ki67_advanced;
        cd.functions.updatePhenotype = null;
        cd.functions.updateVolume = null;
        int Qindex = StandardModels.Ki67_advanced.findPhaseIndex( PhysiCellConstants.Ki67_negative );

        List<double[]> cellPositions = createSphere( sampleСellRadius / 2, sphereRadius );
        for( int i = 0; i < cellPositions.size(); i++ )
        {
            Cell pCell = Cell.createCell( cd, m, VectorUtil.newSum( tumorCenter, cellPositions.get( i ) ) );
            // pCell->functions.volume_update_function=empty_function;
            // pCell->functions.update_phenotype=do_nothing;
            pCell.phenotype.cycle.data.currentPhaseIndex = Qindex;
            pCell.setTotalVolume( volume );
        }
        //        for( int i = 250; i < 260; i += 1 )
        //            visualizer.draw( microenvironment, i, i, resultPath + "/slice_" + i + ".png" );
        for( Visualizer listener : resultListeners )
            listener.init();

        double startSimulation = System.currentTimeMillis();
        double t = 0.0;
        double nextOutputTime = 0;
        while( t < tMax )
        {
            if( Math.abs( t - nextOutputTime ) < 0.001 )
            {
                if( outputReport )
                    writeCellReport( m, resultPath + "/Report_" + (int)t + ".txt" );
                for( Visualizer listener : resultListeners )
                    listener.saveResult( m, t );
                nextOutputTime += outputInterval;
            }
            ( (CellContainer)m.agentContainer ).updateAllCells( m, t, dt, dt, dt );
            t += dt;
        }
        if( outputReport )
            writeCellReport( m, resultPath + "/Report_" + (int)t + ".txt" );
        for( Visualizer listener : resultListeners )
            listener.finish();
        System.out.println( "Done. Elapsed time: " + ( System.currentTimeMillis() - startSimulation ) / 1000 );
    }

    static void writeCellReport(Microenvironment m, String fileName) throws Exception
    {
        Set<BasicAgent> agents = m.getAgents();
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( fileName ) ) ))
        {
            for( BasicAgent agent : agents )
            {
                Cell cell = (Cell)agent;
                String phenotype = cell.phenotype.cycle.currentPhase().name;
                bw.append( cell.ID + "\t" + format( cell.position[0] ) + "\t" + format( cell.position[1] ) + "\t"
                        + format( cell.position[2] ) + "\t" );
                bw.append( format( cell.phenotype.geometry.radius ) + "\t" + format( cell.phenotype.volume.total ) + "\t"
                        + format( cell.phenotype.volume.nuclear_fluid ) + "\t" + format( cell.phenotype.volume.nuclear_solid ) + "\t"
                        + format( cell.phenotype.volume.cytoplasmic_fluid ) + "\t" + format( cell.phenotype.volume.cytoplasmic_solid )
                        + "\t" + format( cell.phenotype.volume.calcified_fraction ) + "\t" + phenotype +
                        // "\t"+ cell.phenotype.cycle.phases[cell.phenotype.current_phase_index].elapsed_time +std::endl;       
                        "\t" + cell.phenotype.cycle.data.elapsedTimePhase + "\n" );
            }
        }
    }

    public static List<double[]> createSphere(double cellRadius, double sphereRadius)
    {
        List<double[]> cells = new ArrayList<>();
        int xc = 0, yc = 0, zc = 0;
        double x_spacing = cellRadius * Math.sqrt( 3 );
        double y_spacing = cellRadius * 2;
        double z_spacing = cellRadius * Math.sqrt( 3 );


        // std::vector<double> cylinder_center(3,0.0);

        for( double z = -sphereRadius; z < sphereRadius; z += z_spacing, zc++ )
            for( double x = -sphereRadius; x < sphereRadius; x += x_spacing, xc++ )
                for( double y = -sphereRadius; y < sphereRadius; y += y_spacing, yc++ )
                {
                    double[] tempPoint = new double[3];
                    tempPoint[0] = x + ( zc % 2 ) * 0.5 * cellRadius;
                    tempPoint[1] = y + ( xc % 2 ) * cellRadius;
                    tempPoint[2] = z;

                    if( Math.sqrt( VectorUtil.norm_squared( tempPoint ) ) < sphereRadius )
                        cells.add( tempPoint );
                }
        return cells;
    }

    public static String format(double val)
    {
        return format.format( val );
    }
}