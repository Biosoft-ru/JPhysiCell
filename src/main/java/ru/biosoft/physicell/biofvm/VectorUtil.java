package ru.biosoft.physicell.biofvm;

import java.util.Arrays;

import ru.biosoft.physicell.core.PhysiCellUtilities;

/*
#############################################################################
# If you use BioFVM in your project, please cite BioFVM and the version     #
# number, such as below:                                                    #
#                                                                           #
# We solved the diffusion equations using BioFVM (Version 1.1.7) [1]        #
#                                                                           #
# [1] A. Ghaffarizadeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730 #
#                                                                           #
#############################################################################
#                                                                           #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)   #
#                                                                           #
# Copyright (c) 2015-2017, Paul Macklin and the BioFVM Project              #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions are    #
# met:                                                                      #
#                                                                           #
# 1. Redistributions of source code must retain the above copyright notice, #
# this list of conditions and the following disclaimer.                     #
#                                                                           #
# 2. Redistributions in binary form must reproduce the above copyright      #
# notice, this list of conditions and the following disclaimer in the       #
# documentation and/or other materials provided with the distribution.      #
#                                                                           #
# 3. Neither the name of the copyright holder nor the names of its          #
# contributors may be used to endorse or promote products derived from this #
# software without specific prior written permission.                       #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED #
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           #
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER #
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       #
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#                                                                           #
#############################################################################
*/
public class VectorUtil
{
    /**
    * Creates new array with values v1[i] - v2[i]
    */
    public static double[] newDiff(double[] v1, double[] v2)
    {
        double[] v = new double[v1.length];
        for( int i = 0; i < v1.length; i++ )
            v[i] = v1[i] - v2[i];
        return v;
    }

    /**
     * Creates new array with values v1[i] + v2[i]
     */
    public static double[] newSum(double[] v1, double[] v2)
    {
        double[] v = new double[v1.length];
        for( int i = 0; i < v1.length; i++ )
            v[i] = v1[i] + v2[i];
        return v;
    }

    /**
     * Creates new array with values v1[i] * v2[i]
     */
    public static double[] newProd(double[] v1, double[] v2)
    {
        double[] v = new double[v1.length];
        for( int i = 0; i < v1.length; i++ )
            v[i] = v1[i] * v2[i];
        return v;
    }

    /**
     * Creates new array with values v1[i] * d
     */
    public static double[] newProd(double[] v1, double d)
    {
        double[] v = new double[v1.length];
        for( int i = 0; i < v1.length; i++ )
            v[i] = v1[i] * d;
        return v;
    }

    /**
     * Creates new array with values v1[i] / v2[i]
     */
    public static double[] newDiv(double[] v1, double[] v2)
    {
        double[] v = new double[v1.length];
        for( int i = 0; i < v1.length; i++ )
            v[i] = v1[i] / v2[i];
        return v;
    }

    /**
     * Returns new array with values v1[i] * d
     */
    public static double[] newProd(double d, double[] v1)
    {
        double[] v = new double[v1.length];
        for( int i = 0; i < v1.length; i++ )
            v[i] = v1[i] * d;
        return v;
    }

    /**
     * Returns new array which contains each element of v1 increased by value d
     */
    public static double[] newSum(double d, double[] v1)
    {
        double[] v = new double[v1.length];
        for( int i = 0; i < v1.length; i++ )
            v[i] = v1[i] + d;
        return v;
    }

    /**
    * Returns new array which contains each element of v1 increased by value d
    */
    public static double[] newSum(double[] v1, double d)
    {
        double[] v = new double[v1.length];
        for( int i = 0; i < v1.length; i++ )
            v[i] = v1[i] + d;
        return v;
    }

    /**
    * Modifies array v1: new ith value is d - v1[i]
    */
    public static void diff(double d, double[] v)
    {
        for( int i = 0; i < v.length; i++ )
            v[i] = d - v[i];
    }

    /**
     * Creates new array with values v1[i] - d
     */
    public static double[] newDiff(double[] v1, double d)
    {
        double[] v = new double[v1.length];
        for( int i = 0; i < v1.length; i++ )
            v[i] = v1[i] - d;
        return v;
    }

    /**
    * Modifies array v1: new ith value is v1[i]+v2[i]
    */
    public static void sum(double[] v1, double[] v2)
    {
        for( int i = 0; i < v1.length; i++ )
            v1[i] += v2[i];
    }

    /**
     * Modifies v1 array with values v1[i] - v2[i]
     */
    public static void diff(double[] v1, double[] v2)
    {
        for( int i = 0; i < v1.length; i++ )
            v1[i] -= v2[i];
    }

    /**
    * Modifies v1 array with values v1[i] / v2[i]
    */
    public static void div(double[] v1, double[] v2)
    {
        for( int i = 0; i < v1.length; i++ )
            v1[i] /= v2[i];
    }

    /**
     * Modifies v1 array with values v1[i] * a
     */
    public static void prod(double[] v1, double a)
    {
        for( int i = 0; i < v1.length; i++ )
            v1[i] *= a;
    }

    /**
     * Modifies v1 array with values v1[i] * v2[i]
     */
    public static void prod(double[] v1, double[] v2)
    {
        for( int i = 0; i < v1.length; i++ )
            v1[i] *= v2[i];
    }

    /**
     * Modifies v1 array with values v1[i] / a
     */
    public static void div(double[] v1, double a)
    {
        for( int i = 0; i < v1.length; i++ )
            v1[i] /= a;
    }

    /**
     * Returns String representing given array: "v[0] v[1] ... v[n]"
     */
    public static String print(double[] v)
    {
        StringBuilder sb = new StringBuilder();
        for( int i = 0; i < v.length; i++ )
        {
            sb.append( v[i] );
            sb.append( " " );
        }
        return sb.toString();
    }

    /**
     * Returns String representing given array: "v[0] v[1] ... v[n]"
     */
    public static String print(int[] v)
    {
        StringBuilder sb = new StringBuilder();
        for( int i = 0; i < v.length; i++ )
        {
            sb.append( v[i] );
            sb.append( " " );
        }
        return sb.toString();
    }

    /**
     * Returns String representing given array: "v[0] v[1] ... v[n]"
     */
    public static String print(String[] v)
    {
        StringBuilder sb = new StringBuilder();
        for( int i = 0; i < v.length; i++ )
        {
            sb.append( v[i] );
            sb.append( " " );
        }
        return sb.toString();
    }

    /**
     * Returns new normalized vector i.e. divided by vector norm: [ v[0]/norm v[1]/norm ... v[n]/norm)<br>
     * Where norm = sqrt( v[0]^2 + v[1]^2 + ... + v[n]^2 ) 
     */
    public static double[] newNormalize(double[] v)
    {
        double norm = 0.0;
        for( int i = 0; i < v.length; i++ )
        {
            norm += v[i] * v[i];
        }
        norm = Math.sqrt( norm );

        if( norm < 1E-16 )
            return new double[v.length];
        // If the norm is small, normalizing doens't make sense. 
        // Just set the entire vector to zero. 
        //        boolean iWarnedYou = false;
        //        if( norm < 1E-16 )
        //        {
        //            if( !iWarnedYou )
        //            {
        //                System.out.println( "Warning and FYI: Very small vector are normalized to 0 vector" );
        //                iWarnedYou = true;
        //            }
        //            return new double[v.length];
        //        }

        double[] output = new double[v.length];
        for( int i = 0; i < v.length; i++ )
        {
            output[i] = v[i] / norm;
        }
        return output;
    }

    /**
     * Normalizes given vector i.e. divide each element by vector norm: [ v[0]/norm v[1]/norm ... v[n]/norm)<br>
     * Where norm = sqrt( v[0]^2 + v[1]^2 + ... + v[n]^2 ) 
     */
    public static void normalize(double[] v)
    {
        double norm = 1e-32;

        for( int i = 0; i < v.length; i++ )
        {
            norm += v[i] * v[i];
        }
        norm = Math.sqrt( norm );
        for( int i = 0; i < v.length; i++ )
        {
            v[i] /= norm;
        }
        // If the norm is small, normalizing doens't make sense. 
        // Just set the entire vector to zero.
        if( norm <= 1e-16 )
        {
            for( int i = 0; i < v.length; i++ )
            {
                v[i] = 0;
            }
        }
    }

    /**
     * Returns squared norm of the vector,<br>
     * Where norm = v[0]^2 + v[1]^2 + ... + v[n]^2 
     */
    public static double norm_squared(double[] v)
    {
        double out = 0;
        for( int i = 0; i < v.length; i++ )
            out += v[i] * v[i];
        return out;
    }

    /**
     * Returns norm of the vector,<br>
     * Where norm = sqrt( v[0]^2 + v[1]^2 + ... + v[n]^2 ) 
     */
    public static double norm(double[] v)
    {
        return Math.sqrt( norm_squared( v ) );
    }

    /**
     * Returns distance between two vector,<br>
     */
    public static double dist(double[] p1, double[] p2)
    {
        return VectorUtil.norm( VectorUtil.newDiff( p1, p2 ) );
    }

    /**
     * Returns greatest in modulus value of the vector:<br>
     * result = v[k]: |v[k]| >= v[i] for all i
     */
    public double maxabs(double[] v)
    {
        double out = 0.0;
        for( int i = 0; i < v.length; i++ )
        {
            if( Math.abs( v[i] ) > out )
            {
                out = v[i];
            }
        }
        return out;
    }


    /**
     * Returns greatest in modulus value of the diffrences between two arrays components:<br>
     * result = v1[k]-v2[k]: |result| >= v1[i] - v2[i] for all i
     */
    public static double max_abs_difference(double[] v1, double[] v2)
    {
        double out = 0.0;
        for( int i = 0; i < v1.length; i++ )
        {
            if( Math.abs( v1[i] - v2[i] ) > out )
            {
                out = v1[i] - v2[i];
            }
        }
        return out;
    }

    /**
     * Returns new array which is exponentiated given array
     * result[i] = EXP( exponent[i] )
     */
    public static double[] exponentiate(double[] exponent)
    {
        double[] out = new double[exponent.length];
        for( int i = 0; i < out.length; i++ )
        {
            out[i] = Math.exp( exponent[i] );
        }

        return out;
    }

    /**
     * Modifies given array: sets all its values to random uniformly distributed numbers between -1 and 1
     */
    public static void randomize(double[] v)
    {
        for( int i = 0; i < v.length; i++ )
        {
            v[i] = -1 + 2 * PhysiCellUtilities.UniformRandom();
        }
    }

    /**
     * Creates and returns random vector of length size, each value is a random number uniformly distributed between min and max
     */
    public static double[] random(int size, double min, double max)
    {
        double[] result = new double[size];
        for( int i = 0; i < size; i++ )
            result[i] = min + ( max - min ) * PhysiCellUtilities.UniformRandom();
        return result;
    }

    /**
     * Modifies y array so that y[i] = y[i] + a*x[i]  
     */
    public static void axpy(double[] y, double a, double[] x)
    {
        for( int i = 0; i < y.length; i++ )
        {
            y[i] += a * x[i];
        }
    }

    /**
     * Modifies y array so that y[i] = y[i] + a[i]*x[i]  
     */
    public static void axpy(double[] y, double[] a, double[] x)
    {
        for( int i = 0; i < y.length; i++ )
        {
            y[i] += a[i] * x[i];
        }
    }

    /**
     * Modifies y array so that y[i] = y[i] - a*x[i]  
     */
    public static void naxpy(double[] y, double a, double[] x)
    {
        for( int i = 0; i < y.length; i++ )
        {
            y[i] -= a * x[i];
        }
    }

    /**
     * Modifies y array so that y[i] = y[i] - a[i]*x[i]  
     */
    public static void naxpy(double[] y, double[] a, double[] x)
    {
        for( int i = 0; i < y.length; i++ )
        {
            y[i] -= a[i] * x[i];
        }
    }

    //
    //// turn a delimited character array (e.g., csv) into a vector of doubles
    //
    //void csv_to_vector( const char* buffer , std::vector<double>& vect )
    //{
    //    vect.resize(0); 
    //    unsigned int i=0;
    //    while( i < strlen( buffer )  )
    //    {
    //        // churn through delimiters, whitespace, etc. to reach the next numeric term
    //        while( isdigit( buffer[i] ) == false && buffer[i] != '.' && buffer[i] != '-' && buffer[i] != 'e' && buffer[i] != 'E' )
    //        { i++; } 
    //        char* pEnd; 
    //        if( i < strlen(buffer) ) // add this extra check in case of a final character, e.g., ']'
    //        {
    //            vect.push_back( strtod( buffer+i , &pEnd ) ); 
    //            i = pEnd - buffer; 
    //        }
    //    }           
    //    return; 
    //}
    //    void csv_to_vector( const char* buffer , std::vector<double>& vect )
    //    {
    //        vect.resize(0); 
    //        unsigned int i=0;
    //        while( i < strlen( buffer )  )
    //        {
    //            // churn through delimiters, whitespace, etc. to reach the next numeric term
    //            while( isdigit( buffer[i] ) == false && buffer[i] != '.' && buffer[i] != '-' && buffer[i] != 'e' && buffer[i] != 'E' )
    //            { i++; } 
    //            char* pEnd; 
    //            if( i < strlen(buffer) ) // add this extra check in case of a final character, e.g., ']'
    //            {
    //                vect.push_back( strtod( buffer+i , &pEnd ) ); 
    //                i = pEnd - buffer; 
    //            }
    //        }           
    //        return; 
    //    }

    //
    //char* vector_to_csv( const std::vector<double>& vect )
    //{ 
    //    static int datum_size = 16;  // format = %.7e, 1 (sign) + 1 (lead) + 1 (decimal) + 7 (figs) + 2 (e, sign) + 3 (exponent) + 1 (delimiter) = 16
    //    // this is approximately the same at matlab long for single precision. 
    //    // If you want better precision, use a binary data format like matlab, or (in the future) HDF 
    //
    //    char* buffer; 
    //    buffer = new char[ datum_size * vect.size() ];
    //    
    //    int position = 0; 
    //    for( unsigned int j=0; j < vect.size()-1 ; j++ )
    //    {
    //        position += sprintf( buffer+position , "%.7e," , vect[j] ); 
    //    }
    //    sprintf( buffer + position , "%.7e" , vect[ vect.size()-1 ] ); 
    //    
    //    return buffer; 
    //}
    //
    //void vector_to_csv_safe( const std::vector<double>& vect , char*& buffer )
    //{ 
    //    static int datum_size = 16;  // format = %.7e, 1 (sign) + 1 (lead) + 1 (decimal) + 7 (figs) + 2 (e, sign) + 3 (exponent) + 1 (delimiter) = 16
    //    // this is approximately the same at matlab long for single precision. 
    //    // If you want better precision, use a binary data format like matlab, or (in the future) HDF 
    //
    //    if( buffer )
    //    { delete [] buffer; } 
    //    buffer = new char[ datum_size * vect.size() ];
    //    std::cout << __LINE__ << std::endl; 
    //    
    //    int position = 0; 
    //    for( unsigned int j=0; j < vect.size()-1 ; j++ )
    //    {
    //        position += sprintf( buffer+position , "%.7e," , vect[j] ); 
    //    }
    //    sprintf( buffer + position , "%.7e" , vect[ vect.size()-1 ] ); 
    //    return; 
    //}
    //
    //void vector_to_csv( const std::vector<double>& vect , char*& buffer )
    //{ 
    //    // %.7e is approximately the same at matlab longe for single precision. 
    //    // If you want better precision, use a binary data format like matlab, or (in the future) HDF 
    //
    //    int position = 0; 
    //    for( unsigned int j=0; j < vect.size()-1 ; j++ )
    //    {
    //        position += sprintf( buffer+position , "%.7e," , vect[j] ); 
    //    }
    //    sprintf( buffer + position , "%.7e" , vect[ vect.size()-1 ] ); 
    //    return; 
    //}
    //
    //// turn a delimited character array (e.g., csv) into a vector of doubles
    //
    //void list_to_vector( const char* buffer , std::vector<double>& vect , char delim ) //
    //{
    //    csv_to_vector( buffer , vect ); 
    //    return; 
    //}
    //
    //char* vector_to_list( const std::vector<double>& vect , char delim )
    //{ 
    //    static int datum_size = 16;  // format = %.7e, 1 (sign) + 1 (lead) + 1 (decimal) + 7 (figs) + 2 (e, sign) + 3 (exponent) + 1 (delimiter) = 16
    //    // this is approximately the same at matlab long for single precision. 
    //    // If you want better precision, use a binary data format like matlab, or (in the future) HDF 
    //
    //    char* buffer; 
    //    buffer = new char[ datum_size * vect.size() ];
    //    
    //    int position = 0; 
    //    for( unsigned int j=0; j < vect.size()-1 ; j++ )
    //    {
    //        position += sprintf( buffer+position , "%.7e%c" , vect[j] , delim ); 
    //    }
    //    sprintf( buffer + position , "%.7e" , vect[ vect.size()-1 ] ); 
    //    
    //    return buffer; 
    //}
    //
    //void vector_to_list_safe( const std::vector<double>& vect , char*& buffer , char delim )
    //{ 
    //    static int datum_size = 16;  // format = %.7e, 1 (sign) + 1 (lead) + 1 (decimal) + 7 (figs) + 2 (e, sign) + 3 (exponent) + 1 (delimiter) = 16
    //    // this is approximately the same at matlab long for single precision. 
    //    // If you want better precision, use a binary data format like matlab, or (in the future) HDF 
    //
    //    if( buffer )
    //    { delete [] buffer; } 
    //    buffer = new char[ datum_size * vect.size() ];
    //    
    //    int position = 0; 
    //    for( unsigned int j=0; j < vect.size()-1 ; j++ )
    //    {
    //        position += sprintf( buffer+position , "%.7e%c" , vect[j] , delim ); 
    //    }
    //    sprintf( buffer + position , "%.7e" , vect[ vect.size()-1 ] ); 
    //    return; 
    //}
    //
    //void vector_to_list( const std::vector<double>& vect , char*& buffer , char delim )
    //{ 
    //    // %.7e is approximately the same at matlab longe for single precision. 
    //    // If you want better precision, use a binary data format like matlab, or (in the future) HDF 
    //
    //    int position = 0; 
    //    for( unsigned int j=0; j < vect.size()-1 ; j++ )
    //    {
    //        position += sprintf( buffer+position , "%.7e%c" , vect[j] , delim ); 
    //    }
    //    sprintf( buffer + position , "%.7e" , vect[ vect.size()-1 ] ); 
    //    return; 
    //}
    //
    //// faster version if you know there are only 3 components
    //void vector3_to_list( const std::vector<double>& vect , char*& buffer , char delim )
    //{ 
    //    // %.7e is approximately the same at matlab longe for single precision. 
    //    // If you want better precision, use a binary data format like matlab, or (in the future) HDF 
    //    sprintf( buffer, "%.7e%c%.7e%c%.7e", vect[0] , delim, vect[1] , delim , vect[2] );
    //    return; 
    //}
    //

    /**
     * Returns dot product of two arrays: a[0]*b[0]+...+a[n]*b[n]
     */
    double dot_product(double[] a, double[] b)
    {
        double out = 0.0;
        for( int i = 0; i < a.length; i++ )
        {
            out += ( a[i] * b[i] );
        }
        return out;
    }

    double[] cross_product(double[] a, double[] b)
    {
        double[] out = new double[3];
        out[0] = a[1] * b[2] - a[2] * b[1];
        out[1] = a[2] * b[0] - a[0] * b[2];
        out[2] = a[0] * b[1] - a[1] * b[0];
        return out;
    }

    /**
     * Returns new array of given size on the base of given array. 
     * If new size is smaller than v.length then new array is simply part of given array
     * Else new array contains v entirely at the start
     * If initial array is null - return array of length 1
     */
    public static boolean[] resize(boolean[] v, int size)
    {
        if( v == null )
            return new boolean[size];
        return Arrays.copyOf( v, size );
    }

    /**
     * Returns new array of given size on the base of given array.<br>
     * If new size is smaller than v.length then new array is simply part of given array<br>
     * Else new array contains v entirely at the start<br>
     * If initial array is null - return array of length 1
     */
    public static double[] resize(double[] v, int size)
    {
        if( v == null )
            return new double[size];
        return Arrays.copyOf( v, size );
    }

    /**
     * Returns new array of given size on the base of given array.<br>
     * If new size is smaller than v.length then new array is simply part of given array<br>
     * Else new array contains v entirely at the start, all other values are equal to val<br>
     * If initial array is null - return array of length 1 with result[0] = val
     */
    public static double[] resize(double[] v, int size, double val)
    {
        double[] result = new double[size];
        if( v == null || v.length == 0 )
        {
            for( int i = 0; i < size; i++ )
                result[i] = val;
            return result;
        }
        if( size <= v.length )
            return Arrays.copyOf( v, size );
        System.arraycopy( v, 0, result, 0, v.length );
        for( int i = v.length; i < size; i++ )
            result[i] = val;
        return result;

    }

    /**
     * Returns new array of given size on the base of given array.<br>
     * If new size is smaller than v.length then new array is simply part of given array<br>
     * Else new array contains v entirely at the start, all other elements are arrays of the same size as v[0]
     */
    public static int[][] resize(int[][] v, int size)
    {
        if( v == null || v.length == 0 )
            return new int[size][];
        if( size <= v.length )
            return Arrays.copyOf( v, size );
        int[][] result = Arrays.copyOf( v, size );
        int innerSize = v[0] == null ? 0 : v[0].length;
        for( int i = v.length; i < result.length; i++ )
            result[i] = new int[innerSize];
        return result;
    }

    /**
     * Returns new array of given size on the base of given array.<br>
     * If new size is smaller than v.length then new array is simply part of given array<br>
     * Else new array contains v entirely at the start, all other elements are arrays of the same size as v[0]
     */
    public static double[][] resize(double[][] v, int size)
    {
        if( v == null || v.length == 0 )
            return new double[size][];
        if( size <= v.length )
            return Arrays.copyOf( v, size );
        double[][] result = Arrays.copyOf( v, size );
        int innerSize = v[0].length;
        for( int i = v.length; i < result.length; i++ )
            result[i] = new double[innerSize];
        return result;
    }

    /**
     * Returns new array which length is arr.length +1. Such that:<br>
     * for i < arr.length result[i] = arr[i]<br>
     * result[arr.length] = val
     */
    public static <T> T[] push_back(T[] arr, T val)
    {
        T[] result = Arrays.copyOf( arr, arr.length + 1 );
        result[result.length - 1] = val;
        return result;
    }

    /**
     * Returns new array which length is arr.length +1. Such that:<br>
     * for i < arr.length result[i] = arr[i]<br>
     * result[arr.length] = val<br>
     * If initial array is null or empty returns array with one value: val
     */
    public static boolean[] push_back(boolean[] arr, boolean val)
    {
        if( arr == null || arr.length == 0 )
            return new boolean[] {val};
        boolean[] result = Arrays.copyOf( arr, arr.length + 1 );
        result[result.length - 1] = val;
        return result;
    }

    /**
     * Returns new array which length is arr.length +1. Such that:<br>
     * for i < arr.length result[i] = arr[i]<br>
     * result[arr.length] = val<br>
     * If initial array is null or empty returns array with one value: val
     */
    public static double[] push_back(double[] arr, double val)
    {
        if( arr == null || arr.length == 0 )
            return new double[] {val};
        double[] result = Arrays.copyOf( arr, arr.length + 1 );
        result[result.length - 1] = val;
        return result;
    }

    /**
     * Returns new array which length is arr.length +1. Such that:<br>
     * for i < arr.length result[i] = arr[i]<br>
     * result[arr.length] = val<br>
     * If initial array is null or empty returns array with one value: val
     */
    public static int[] push_back(int[] arr, int val)
    {
        if( arr == null || arr.length == 0 )
            return new int[] {val};
        int[] result = Arrays.copyOf( arr, arr.length + 1 );
        result[result.length - 1] = val;
        return result;
    }

    /**
     * Returns new array of size length with all values equal to template:<br>
     * result[i] = template
     */
    public static double[][] assign(int length, double[] template)
    {
        double[][] result = new double[length][];
        for( int i = 0; i < length; i++ )
            result[i] = Arrays.copyOf( template, template.length );
        return result;
    }

    /**
     * Returns new array of size length with all values equal to template:<br>
     * result[i] = template
     */
    public static boolean[][] assign(int length, boolean[] template)
    {
        boolean[][] result = new boolean[length][];
        for( int i = 0; i < length; i++ )
            result[i] = Arrays.copyOf( template, template.length );
        return result;
    }

    /**
     * Returns new array of size length with all values equal to template:<br>
     * result[i] = template
     */
    public static double[] assign(int length, double template)
    {
        double[] result = new double[length];
        for( int i = 0; i < length; i++ )
            result[i] = template;
        return result;
    }

    /**
     * Returns true if two arrays are of the same size and all corresponding elements are equal: arr1[i] = arr2[i]
     */
    public static boolean equals(double[] arr1, double[] arr2)
    {
        if( arr1 == null || arr2 == null || arr1.length != arr2.length )
            return false;
        for( int i = 0; i < arr1.length; i++ )
        {
            if( arr1[i] != arr2[i] )
                return false;
        }
        return true;
    }

    /**
     * turn a delimited character array (e.g., csv) into a vector of doubles
     */
    public static double[] csv_to_vector(String buffer)
    {
        String[] splited = buffer.split( " " );
        double[] result = new double[splited.length];
        for( int i = 0; i < splited.length; i++ )
            result[i] = Double.parseDouble( splited[i] );
        return result;
        //        List<Double> result = new ArrayList<>();
        //        int length = buffer.length();
        //        for( int i = 0; i < length; i++ )
        //        {
        //            char c = buffer.charAt( i );
        //            // churn through delimiters, whitespace, etc. to reach the next numeric term
        //            if( !Character.isDigit( c ) && c != '.' && c != '-' && c != 'e' && c != 'E' )
        //            {
        //                continue;
        //            }
        //            String pEnd;
        //            if( i < length ) // add this extra check in case of a final character, e.g., ']'
        //            {
        //                result.add( Double.parseDouble( c ) ); //vect.push_back( strtod( buffer+i , &pEnd ) ); 
        //                i = pEnd - buffer; 
        //            }
        //        }  
        //        return result.stream().mapToDouble(d -> d).toArray();
    }
}