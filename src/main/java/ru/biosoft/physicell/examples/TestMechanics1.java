package ru.biosoft.physicell.examples;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellContainer;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.StandardModels;
import ru.biosoft.physicell.ui.Visualizer;
import ru.biosoft.physicell.ui.Visualizer.Section;

public class TestMechanics1
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

    public static double dt = 0.05;
    public static double t_output_interval = 0.05; // 1.0; 
    public static double t_max = 5;
    private static int size = 500;
    private static String resultPath = "C:/Users/Damag/BIOFVM/Mechanics1/";
    private static String resultName = "result";
    private static double[] point1 = new double[] {100, 100, 100};
    private static double[] point2 = new double[] {105.73, 100, 100};
    private static double volume = 4188.790204786391;

    private static Visualizer visualizer = new Visualizer( resultPath, resultName, Section.Z, 100 );
    static
    {
        visualizer.setSaveImage( false );
        visualizer.setDrawDensity( false );
    }

    public static void main(String ... args) throws Exception
    {
        Microenvironment m = new Microenvironment( "substrate scale", size, 20, "minutes", "microns" );
        m.setDensity( 0, "oxygen", "mmHg", 0, 0 );
        CellContainer.createCellContainer( m, 30 );
        CellDefinition cd = StandardModels.createFromDefault( "tumor cell", 0, m );
        CellDefinition.registerCellDefinition( cd );
        cd.phenotype.cycle = StandardModels.Ki67_advanced;
        cd.functions.updatePhenotype = null;
        cd.functions.updateVolume = null;

        int Q_index = StandardModels.Ki67_advanced.findPhaseIndex( PhysiCellConstants.Ki67_negative );

        Cell pCell1 = Cell.createCell( cd, m, point1 );
        pCell1.phenotype.cycle.data.currentPhaseIndex = Q_index;
        /* NOTE: for this experiment, you need to disable volume update function 
         to make sure that volume change are not affecting the distance we measure for the cells.*/
        pCell1.functions.updateVolume = null;
        pCell1.setTotalVolume( volume );

        Cell pCell2 = Cell.createCell( cd, m, point2 );
        pCell2.phenotype.cycle.data.currentPhaseIndex = Q_index;
        pCell2.functions.updateVolume = null;
        pCell2.setTotalVolume( volume );

        pCell1.functions.updateVelocity.execute( pCell1, pCell1.phenotype, dt );
        pCell2.functions.updateVelocity.execute( pCell2, pCell2.phenotype, dt );

        pCell1.setPreviousVelocity( pCell1.velocity[0], pCell1.velocity[1], pCell1.velocity[2] );
        pCell2.setPreviousVelocity( pCell2.velocity[0], pCell2.velocity[1], pCell2.velocity[2] );

        //        for( int i = 0; i < 10; i++ )
        //        {
        //            VectorUtil.axpy( pCell1.position, ( dt / 10.0 ), pCell1.velocity );
        //            VectorUtil.axpy( pCell2.position, ( dt / 10.0 ), pCell2.velocity );
        //            t += dt / 10.0;
        //        }

        double t = 0.0;
        double t_next_output_time = 0;
        visualizer.init();
        while( t < t_max )
        {
            if( Math.abs( t - t_next_output_time ) < dt / 10.0 )
            {
                System.out.println( t + "\t" + VectorUtil.dist( pCell1.position, pCell2.position ) );
                visualizer.saveResult( m, t );
                t_next_output_time += t_output_interval;
            }
            ( (CellContainer)m.agentContainer ).updateAllCells( m, t, dt, dt, dt );
            t += dt;
        }
        visualizer.finish();
    }
}