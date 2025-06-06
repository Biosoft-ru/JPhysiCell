package ru.biosoft.physicell.sample_projects.pred_prey_farmer;

import java.io.InputStream;

import ru.biosoft.physicell.biofvm.ConstantCoefficientsLOD3D;
import ru.biosoft.physicell.core.CellContainer;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.ui.GIFGenerator;
import ru.biosoft.physicell.ui.render.Visualizer3D;
import ru.biosoft.physicell.xml.ModelReader;

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
public class Main
{
    private static String settingsPath = "config/PhysiCell_settings.xml";
    private static String resultPath = "C:/Users/Damag/BIOFVM/projects/perd_prey_farmer/October";

    public static void main(String ... strings) throws Exception
    {
        if( strings != null && strings.length > 0 )
            resultPath = strings[0];

        InputStream settings = Main.class.getResourceAsStream( settingsPath );
        Model model = new ModelReader().read( settings, PredPreyFarmer.class );
        double mechanics_voxel_size = 30;
        ConstantCoefficientsLOD3D solver = new ConstantCoefficientsLOD3D();
        solver.setPrallel( false );
        model.getMicroenvironment().setSolver( solver );
        model.setResultFolder( resultPath );
        model.createContainer( mechanics_voxel_size, CellContainer.DEFAULT_NAME );
        model.setWriteDensity( false );
        model.setSaveFull( true );
        model.setSaveImg( true );

        Visualizer3D visualizer = new Visualizer3D( resultPath, "3d", model.getMicroenvironment() );
        visualizer.addResultGenerator( new GIFGenerator( resultPath, "3d.gif" ) );
        model.addVisualizer( visualizer );
//        model.addGIFVisualizer( 0, "food" ).setStubstrateIndex( 0 ).setMaxDensity( 10 ).setAgentVisualizer( new PPFVisualizer() );
        //        model.addGIFVisualizer( 0, "prey signal" ).setStubstrateIndex( 1 ).setMaxDensity( 10 ).setAgentVisualizer( new PPFVisualizer() );
        //        model.addGIFVisualizer( 0, "predator signal" ).setStubstrateIndex( 2 ).setMaxDensity( 10 ).setAgentVisualizer( new PPFVisualizer() );
        model.init();
        System.out.println( model.display() );
        model.simulate();
    }
}