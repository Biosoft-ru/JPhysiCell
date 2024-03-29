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
public class Geometry implements Cloneable
{
    public double radius;
    private double nuclearRadius;
    double surfaceArea;
    public double polarity;
    static double four_thirds_pi = 4.188790204786391;
    static double the_constant = 4.835975862049409; // 4pi / (4pi/3)^(2/3)

    public Geometry()
    {
        // reference values for MCF-7, based on volume = 2494 cubic microns nuclear volume = 540 cubic microns 
        radius = 8.412710547954228;
        nuclearRadius = 5.051670902881889;
        surfaceArea = 889.3685284131693;
        polarity = 0.0;
    }

    void updateRadius(Cell pCell, Phenotype phenotype, double dt)
    {
        radius = phenotype.volume.total;
        radius /= four_thirds_pi;
        radius = Math.pow( radius, 0.333333333333333333333333333333333333333 );
    }

    public double getNuclearRadius()
    {
        return nuclearRadius;
    }

    public double getRadius()
    {
        return radius;
    }

    void updateNuclearRadius(Cell pCell, Phenotype phenotype, double dt)
    {
        nuclearRadius = phenotype.volume.nuclear;
        nuclearRadius /= four_thirds_pi;
        nuclearRadius = Math.pow( nuclearRadius, 0.333333333333333333333333333333333333333 );
    }

    void updateSurfaceArea(Cell pCell, Phenotype phenotype, double dt)
    {

        surfaceArea = Math.pow( phenotype.volume.total, 0.666666666666667 );
        surfaceArea /= the_constant;
    }

    public void update(Cell pCell, Phenotype phenotype, double dt)
    {
        updateRadius( pCell, phenotype, dt );
        updateNuclearRadius( pCell, phenotype, dt );
        surfaceArea = 3 * phenotype.volume.total / radius; // surface area = 4*pi*r^2 = (4/3)*pi*r^3 / (r/3) 
        //        surfaceArea /= radius;
        //        surfaceArea *= 3.0;
    }

    @Override
    public Geometry clone()
    {
        try
        {
            return (Geometry)super.clone();
        }
        catch( CloneNotSupportedException e )
        {
            throw new InternalError( e );
        }
    }

}