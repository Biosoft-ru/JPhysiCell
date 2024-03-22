package ru.biosoft.physicell.core;

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
public class CellFunctions
{
    Model model;
    public instantiate_cell instantiate_cell;
    public VolumeUpdate updateVolume;
    public UpdateMigrationBias updateMigration;
    public CustomCellRule customCellRule;
    public UpdatePhenotype updatePhenotype;
    public pre_update_intracellular pre_update_intracellular;
    public post_update_intracellular post_update_intracellular;
    public UpdateVelocity updateVelocity;
    public MembraneInteractions membraneInteraction;
    public DistanceCalculator membraneDistanceCalculator;
    public set_orientation set_orientation;
    public Contact contact;

    public CellFunctions clone()
    {

        CellFunctions result = new CellFunctions();
        try
        {
            result.instantiate_cell = instantiate_cell;// == null ? null : instantiate_cell.getClass().newInstance();
            result.updateVolume = updateVolume;// == null ? null : updateVolume.getClass().newInstance();
            result.updateMigration = updateMigration;// == null ? null : updateMigration.getClass().newInstance();
            result.updatePhenotype = updatePhenotype;// == null ? null : updatePhenotype.clone();
            result.pre_update_intracellular = pre_update_intracellular;// == null ? null : pre_update_intracellular.getClass().newInstance();
            result.post_update_intracellular = post_update_intracellular;// == null ? null;
            //                    : post_update_intracellular.getClass().newInstance();
            result.updateVelocity = updateVelocity;// == null ? null : updateVelocity.getClass().newInstance();
            result.membraneInteraction = membraneInteraction;// == null ? null : membraneInteraction.getClass().newInstance();
            result.membraneDistanceCalculator = membraneDistanceCalculator;// == null ? null
            //                    : membraneDistanceCalculator.getClass().newInstance();
            result.set_orientation = set_orientation;// == null ? null : set_orientation.getClass().newInstance();
            result.contact = contact;// == null ? null : contact.getClass().newInstance();
            result.customCellRule = customCellRule;
        }
        catch( Exception ex )
        {
            ex.printStackTrace();
        }
        return result;
    }

    public static abstract class Function implements Cloneable
    {
        public UpdatePhenotype clone()
        {
            try
            {
                return (UpdatePhenotype)super.clone();
            }
            catch( Exception ex )
            {
                return null;
            }
        }

        public String display()
        {
            return getClass().getSimpleName();
        }
    }

    public static abstract class OneCellFunction extends Function
    {
        public abstract void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception;
    }

    public static abstract class VolumeUpdate extends OneCellFunction
    {
    }

    public static abstract class UpdateMigrationBias extends OneCellFunction
    {
    }

    public static abstract class UpdatePhenotype extends OneCellFunction
    {
    }

    public static abstract class UpdateVelocity extends OneCellFunction
    {
    }

    public static abstract class CustomCellRule extends OneCellFunction
    {

    }

    public static abstract class pre_update_intracellular extends OneCellFunction
    {
    }

    public static abstract class post_update_intracellular extends OneCellFunction
    {
    }

    public static abstract class MembraneInteractions extends OneCellFunction
    {
    }

    public static abstract class DistanceCalculator extends Function
    {
        public abstract double execute(Cell pCell, Phenotype phenotype, double dt) throws Exception;
    }

    public static abstract class set_orientation extends OneCellFunction
    {
    }

    public static abstract class Contact extends Function
    {
        public abstract void execute(Cell pCell, Phenotype phenotype, Cell cell2, Phenotype phenotype2, double dt);
    }

    @FunctionalInterface
    public static interface instantiate_cell
    {
        public Cell execute();
    }

    @Override
    public String toString()
    {
        return display();
    }

    public String display()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "Key functions: \n--------------------------------" );
        if( updateMigration != null )
            sb.append( "\n\tMigration bias rule: " + updateMigration.display() );
        if( customCellRule != null )
            sb.append( "\n\tCustom rule: " + customCellRule.display() );
        if( updatePhenotype != null )
            sb.append( "\n\tPhenotype rule: " + updatePhenotype.display() );
        if( updateVolume != null )
            sb.append( "\n\tVolume update function: " + updateVolume.display() );
        if( updateVelocity != null )
            sb.append( "\n\tMechanics function: " + updateVelocity.display() );
        if( contact != null )
            sb.append( "\n\tContact function: " + contact.display() );
        if( membraneDistanceCalculator != null )
            sb.append( "\n\tMembrane distance: " + membraneDistanceCalculator.display() );
        return sb.toString();
    }
}