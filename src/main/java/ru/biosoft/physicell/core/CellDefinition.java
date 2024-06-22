package ru.biosoft.physicell.core;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.CustomCellData.VectorVariable;

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
# Copyright (c) 2015-2022, Paul Macklin and the PhysiCell Project             #
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
public class CellDefinition
{
    public int type;
    public String name;

    public boolean isMovable;

    private Microenvironment m;
    public CellParameters parameters = new CellParameters();
    public CustomCellData custom_data = new CustomCellData();
    public CellFunctions functions = new CellFunctions();
    public Phenotype phenotype = new Phenotype();

    public CellDefinition()
    {
        isMovable = true;
        //        parameters.pReference_live_phenotype = phenotype; //TODO: check
        // set up the default custom data, the default Custom_Cell_Data constructor should take care of this
        // set up the default functions 
        functions.instantiate_cell = null;
        functions.updateVolume = null; // standard_volume_update_function;
        functions.updateMigration = null;

        functions.updatePhenotype = null;
        functions.customCellRule = null;

        functions.updateVelocity = null; // standard_update_cell_velocity;
        functions.membraneInteraction = null;
        functions.membraneDistanceCalculator = null;
        functions.set_orientation = null;
    }

    public CellDefinition(Microenvironment m, int type, String name)
    {
        this();
        this.m = m;
        this.type = type;
        this.name = name;
        phenotype.sync( m );
    }

    public CellDefinition clone(String name, int type, Microenvironment m)
    {
        CellDefinition result = new CellDefinition( m, type, name );
        result.isMovable = this.isMovable;
        result.parameters = parameters.clone();
        result.custom_data = custom_data.clone();
        result.phenotype = phenotype.clone();
        result.functions = functions.clone();
        result.phenotype.sync( m );
        return result;
    }

    public CellDefinition clone(String name, int type)
    {
        CellDefinition result = new CellDefinition( );
        result.name = name;
        result.type = type;
        result.isMovable = this.isMovable;
        result.parameters = parameters.clone();
        result.custom_data = custom_data.clone();
        result.phenotype = phenotype.clone();
        result.functions = functions.clone();
        return result;
    }

    public static void copy(CellDefinition from, CellDefinition to)
    {
        to.isMovable = from.isMovable;
        to.parameters = from.parameters.clone();
        to.custom_data = from.custom_data.clone();
        to.phenotype = from.phenotype.clone();
        to.functions = from.functions.clone();
    }

    @Override
    public String toString()
    {
        return name + " ( " + type + " ) ";
    }

    public Microenvironment getMicroenvironment()
    {
        return m;
    }

    public String display()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "================================" );
        sb.append( "\n\t" + name + " (" + type + ")" );
        sb.append( "\n================================\n" );

        if( phenotype.cycle != null )
            sb.append( "\n" + phenotype.cycle.display() + "\n" );
        else
        {
            sb.append( "\n\nCycle Model not initialized." );
            sb.append( "\n--------------------------------" );
        }

        sb.append( "\n" + phenotype.death.display() + "\n" );
        sb.append( "\n" + phenotype.motility.display() + "\n" );
        sb.append( "\n" + phenotype.secretion.display() + "\n" );
        sb.append( "\n" + phenotype.cellInteractions.display() + "\n" );
        sb.append( "\n" + phenotype.cellTransformations.display() + "\n" );
        sb.append( "\n" + phenotype.mechanics.display() + "\n" );
        sb.append( "\n" + functions.display() + "\n" );

        //        if( phenotype.intracellular != null )
        //        {
            //                    phenotype.intracellular.display();
        //        }

        CustomCellData pCCD = ( custom_data );
        sb.append( "\nCustom data: " );
        sb.append( "\n--------------------------------" );
        for( Variable var : pCCD.variables )
        {
            if( var.value != 0 )
                sb.append( "\n\t" + var );
        }
        if( !pCCD.vectorVariables.isEmpty() )
        {
            sb.append( "\n\nCustom vector data: " );
            sb.append( "\n--------------------------------" );
            if( !pCCD.vectorVariables.isEmpty() )
            {
                for( VectorVariable var : pCCD.vectorVariables )
                    sb.append( "\n\t" + var );
            }
        }
        return sb.toString();
    }

    public static String display_ptr_as_bool(Object obj)
    {
        return obj == null ? "false" : "true";
    }
}