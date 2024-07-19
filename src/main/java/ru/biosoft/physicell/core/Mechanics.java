package ru.biosoft.physicell.core;

import ru.biosoft.physicell.biofvm.VectorUtil;

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
public class Mechanics implements Cloneable
{
    public double cellCellAdhesionStrength;
    public double cellBMAdhesionStrength;
    public double cellCellRepulsionStrength;
    public double cellBMRepulsionStrength;
    public double[] cellAdhesionAffinities;

    public double relMaxAdhesionDistance; // this is a multiple of the cell (equivalent) radius
    // double maximum_adhesion_distance; // needed? 

    /* for spring attachments */
    public int maxAttachments;
    public double attachmentElasticConstant;

    public double attachmentRate;
    public double detachmentRate;

    /* to be deprecated */
    public double relMaxAttachmentDistance;
    public double relDetachmentDistance;
    public double maxAttachmentRate;

    public Mechanics()
    {
        cellCellAdhesionStrength = 0.4;
        cellBMAdhesionStrength = 4.0;

        cellCellRepulsionStrength = 10.0;
        cellBMRepulsionStrength = 100.0;

        cellAdhesionAffinities = new double[] {1};

        // this is a multiple of the cell (equivalent) radius
        relMaxAdhesionDistance = 1.25;
        // maximum_adhesion_distance = 0.0; 

        /* for spring attachments */
        maxAttachments = 12;
        attachmentElasticConstant = 0.01;

        attachmentRate = 0; // 10.0 prior ot March 2023
        detachmentRate = 0;

        /* to be deprecated */
        relMaxAttachmentDistance = relMaxAdhesionDistance;
        relDetachmentDistance = relMaxAdhesionDistance;

        maxAttachmentRate = 1.0;
    }

    public void initialize(Model model)
    {
        cellAdhesionAffinities = VectorUtil.resize( cellAdhesionAffinities, model.getDefinitionsCount(), 1.0 );
    }

    public double cell_adhesion_affinity(String type_name, Model model)
    {
        int n = model.getCellDefinition( type_name ).type;
        return cellAdhesionAffinities[n];
    }

    void setFullyHeterotypic(Model model)
    {
        cellAdhesionAffinities = VectorUtil.assign( model.getDefinitionsCount(), 1.0 );
    }

    void setFullyHomotypic(Cell pC)
    {
        cellAdhesionAffinities = new double[pC.getModel().getDefinitionsCount()];
        // now find my type and set to 1 
        //  cell_adhesion_affinity( pC->type_name ) = 1.0; 
    }

    // new on July 29, 2018
    // change the ratio without changing the repulsion strength or equilibrium spacing 
    void setRelMaxAdhesionDistance(double newValue)
    {
        // get old equilibrium spacing, based on equilibriation of pairwise adhesive/repulsive forces at that distance. 
        // relative equilibrium spacing (relative to mean cell radius)
        double temp1 = Math.sqrt( cellCellAdhesionStrength / cellCellRepulsionStrength );
        double temp2 = 1.0 - temp1; //  1 - sqrt( alpha_CCA / alpha_CCR );
        double sRelative = 2 * temp2; // 2*( 1 - sqrt( alpha_CCA / alpha_CCR ) ); 

        temp1 /= relMaxAdhesionDistance; // sqrt( alpha_CCA / alpha_CCR)/f;
        temp2 = 1.0 - temp1; // 1 - sqrt( alpha_CCA / alpha_CCR )/f;
        sRelative /= temp2; // 2*( 1 - sqrt( alpha_CCA / alpha_CCR ) ) / ( 1-1/f) ; 

        // now, adjust the relative max adhesion distance 
        relMaxAdhesionDistance = newValue;

        // adjust the adhesive coefficient to preserve the old equilibrium distance
        temp1 = sRelative / 2.0;
        temp2 = 1.0 - temp1; // 1 - s_relative/2.0 
        temp1 = 1.0 - temp1 / relMaxAdhesionDistance; // 1.0 -s_relative/(2*relative_maximum_adhesion_distance); 
        temp2 /= temp1;
        temp2 *= temp2;
        cellCellAdhesionStrength = cellCellRepulsionStrength * temp2;
    }

    // new on July 29, 2018
    // set the cell-cell equilibrium spacing, accomplished by changing the 
    // cell-cell adhesion strength, while leaving the cell-cell repulsion 
    // strength and the maximum adhesion distance unchanged 
    public void setRelEquilibriumDistance(double newValue)
    {
        if( newValue > 2.0 )
        {
            System.out.println( "Warning in Mechanics function :\n" + "\tAttempted to set equilibrium distance exceeding two cell radii.\n"
                    + "\tWe will cap the equilibrium distance at 2.0 cell radii.\n" + "****\n\n" );
            newValue = 2.0;
        }
        // adjust the adhesive coefficient to achieve the new (relative) equilibrium distance
        double temp1 = newValue / 2.0;
        double temp2 = 1.0 - temp1; // 1 - s_relative/2.0
        temp1 = 1 - temp1 / relMaxAdhesionDistance; // 1.0 -s_relative/(2*relative_maximum_adhesion_distance);         
        temp2 /= temp1;
        temp2 *= temp2;
        cellCellAdhesionStrength = cellCellRepulsionStrength * temp2;
    }

    public void setAbsEquilibriumDistance(Phenotype phenotype, double newValue)
    {
        setRelEquilibriumDistance( newValue / phenotype.geometry.radius );
    }

    @Override
    public Mechanics clone()
    {
        try
        {
            Mechanics result = (Mechanics)super.clone();
            result.cellAdhesionAffinities = cellAdhesionAffinities.clone();
            return result;
        }
        catch( CloneNotSupportedException e )
        {
            throw ( new InternalError( e ) );
        }
    }

    public String display()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "Mechanics:" );
        sb.append( "\n--------------------------------" );
        sb.append( "\n\tcell_cell_adhesion_strength: " + cellCellAdhesionStrength + "\n" + "\tcell_cell_repulsion_strength: "
                + cellCellRepulsionStrength + "\n" + "\trel max adhesion dist: " + relMaxAdhesionDistance + "\n"
                + "\tcell_BM_adhesion_strength: " + cellBMAdhesionStrength + "\n" + "\tcell_BM_repulsion_strength: "
                + cellBMRepulsionStrength + "\n" + "\tattachment_elastic_constant: " + attachmentElasticConstant + "\n"
                + "\tattachment_rate: " + attachmentRate + "\n" + "\tdetachment_rate: " + detachmentRate );
        return sb.toString();
    }
}