package ru.biosoft.physicell.sample_projects.biorobots;

import java.io.InputStream;

import ru.biosoft.physicell.biofvm.ConstantCoefficientsLOD3D;
import ru.biosoft.physicell.core.Model;
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
# Copyright (c) 2015-2023, Paul Macklin and the PhysiCell Project             #
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
    private static String resultPath = "C:/Users/Damag/BIOFVM/projects/biorobots/seq";

    public static void main(String ... strings) throws Exception
    {
        InputStream settings = Main.class.getResourceAsStream( settingsPath );
        Model model = new ModelReader().read( settings, Biorobots.class );
        double mechanics_voxel_size = 30;
        ConstantCoefficientsLOD3D solver = new ConstantCoefficientsLOD3D();
        solver.setPrallel( false );
        model.getMicroenvironment().setSolver( solver );
        model.createContainer( mechanics_voxel_size );
        model.setResultFolder( resultPath );
        model.setSaveFull( false );
        model.setSaveImg( false );
        model.addGIFVisualizer( 0, "figure1" ).setStubstrateIndex( 1 ).setMaxDensity( 0.5 ).setAgentVisualizer( new BiorobotsVisualizer( this ) ) );;
        model.init();
        System.out.println( model.display() );

        double tStart = System.nanoTime();
        model.simulate();
        tStart = System.nanoTime() - tStart;

        //        System.out.println( "Total " + tStart / 1E9 );
        //        System.out.println( "Diffusion " + Model.tDiffusion / 1E9 );
        //        System.out.println( "Secretion " + CellContainer.tSecretion / 1E9 );
        //        System.out.println( "Phenotype " + CellContainer.tPhenotype / 1E9 );
        //        System.out.println( "Velocity " + CellContainer.tVelocity / 1E9 );
        //        System.out.println( "Interaction " + CellContainer.tInteraction / 1E9 );
        //        System.out.println( "Contact " + CellContainer.tContact / 1E9 );
        //        System.out.println( "Custom " + CellContainer.tCustom / 1E9 );
        //        System.out.println( "Attachment " + CellContainer.tAttachment / 1E9 );
        //        System.out.println( "Divide " + CellContainer.tDivide / 1E9 );
        //        System.out.println( "Rest " + CellContainer.tRestAll / 1E9 );
        //        System.out.println( "Gradient " + CellContainer.tGradient / 1E9 );
        //        System.out.println( "Total cell " + CellContainer.tTotal / 1E9 );
    }
}
