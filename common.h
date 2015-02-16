/*
 * common.h
 *
 *  Created on: 10 Feb, 2015
 *      Author: marcos
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <math.h>
#include <map>
#include <cstdlib>   //srand, rand
#include <ctime>     //time
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

#define MEDIAN 0
#define PERCENTIL_25 1
#define PERCENTIL_75 2

#define NO_MODE 0             //NO MODE SELECTED
#define MODE_PERMUTATIONS 1   //RUN the application to look for differential expression in genes
#define MODE_SPLITTER 2       //SPLITS INPUT FILE
#define MODE_FDR 3			  //FDR Estimation


template<class TYPE>
/**
 * \brief Transforms number to string
 */
inline string numberToStr(TYPE number)
{
	ostringstream buff;
	buff.precision(15);
	buff << number;
	return buff.str();
}

class SampleGene
{
	public:
		string sampleName;
		double expressionValue;
};

class Gene
{
	public:
			string geneName;
			vector <SampleGene> vGeneValues;
			vector <bool> vGroup;

			unsigned int hits;
			unsigned int nPermutationsDone;
			double pValue;
			double originalDiff;
			double medianCases;
			double originalMedianControl;
			double originalMedianCases;
			double medianControl;
			double foldChange;

			Gene()
			{
				hits = 0;
				nPermutationsDone = 0;
				pValue = -1;
				originalDiff = -1;
				medianCases = 0;
				medianControl = 0;
				originalMedianCases = -1;
				originalMedianControl = -1;
				foldChange = -1;
				//set the random number seed for the rand() operator
				srand ( unsigned ( time (NULL) ) );

			}

			//> Adds a new sample expression value
			inline void addSample(const string & sampleName, bool group, double expressionValue)
			{
				SampleGene newSample;
				newSample.sampleName = sampleName;
				newSample.expressionValue = expressionValue;
				vGeneValues.push_back(newSample);
				vGroup.push_back(group);
			}

			//> Rearranges the elements in a vector
			void random_shuffle()
			{
				std::random_shuffle ( vGroup.begin(), vGroup.end() );
			}

			//> Returns the difference values of the median between controls and cases
			inline double diffMedian()
			{
				vector <double> vControls;
				vector <double> vCases;

				vector<bool>::const_iterator iterGroup = vGroup.begin();

				for (vector<SampleGene>::const_iterator iter = vGeneValues.begin(); iter != vGeneValues.end(); ++iter)
				{
					if((*iterGroup))
					{
						vControls.push_back((*iter).expressionValue);
					}
					else
					{
						vCases.push_back((*iter).expressionValue);
					}
					++ iterGroup;
				}

			    if(vControls.empty() || vCases.empty()) return 0;
			    else
			    {
			        sort(vControls.begin(), vControls.end());
			        sort(vCases.begin(), vCases.end());

			        //Median Controls
			        if(vControls.size() % 2 == 0) medianControl = (vControls[vControls.size()/2 - 1] + vControls[vControls.size()/2]) / 2;
			        else medianControl = vControls[vControls.size()/2];
			        //Median Cases
			        if(vCases.size() % 2 == 0) medianCases = (vCases[vCases.size()/2 - 1] + vCases[vCases.size()/2]) / 2;
			        else medianCases = vCases[vCases.size()/2];

			        if(originalMedianCases < 0 && originalMedianControl < 0)
			        {
			        	originalMedianCases = medianCases;
			        	originalMedianControl = medianControl;
			        	foldChange = originalMedianCases/originalMedianControl;
			        }

			        return fabsf(medianControl - medianCases);
			    }
			}

			//> Returns the difference values of the percentile between controls and cases
			inline double diffPercentile(double percentile)
			{
				//Percentil must be 0.25 0.75
				vector <double> vControls;
				vector <double> vCases;

				vector<bool>::const_iterator iterGroup = vGroup.begin();

				for (vector<SampleGene>::const_iterator iter = vGeneValues.begin(); iter != vGeneValues.end(); ++iter)
				{
					if((*iterGroup))
					{
						vControls.push_back((*iter).expressionValue);
					}
					else
					{
						vCases.push_back((*iter).expressionValue);
					}

					++ iterGroup;
				}

				if(vControls.empty() || vCases.empty()) return 0;
				else
				{
					//Control percentile
					nth_element(vControls.begin(),vControls.begin() + int(vControls.size()*percentile), vControls.end());
					medianControl = *(vControls.begin()+int(vControls.size()*percentile));

					//Case percentile
					nth_element(vCases.begin(),vCases.begin() + int(vCases.size()*percentile), vCases.end());
					medianCases = *(vCases.begin()+int(vCases.size()*percentile));

			        if(originalMedianCases < 0 && originalMedianControl < 0)
			        {
			        	originalMedianCases = medianCases;
			        	originalMedianControl = medianControl;
			        	foldChange = originalMedianCases/originalMedianControl;
			        }

					return fabsf(medianControl - medianCases);
				}
			}

			//https://www.informit.com/guides/content.aspx?g=cplusplus&seqNum=290
};



class GenePvalue
{
 private:
	string gene;
	string diff;
	string caseSample;
	string controlSample;
	string fold_change;
	string hits;
	string nPermutation;
	double pv;
	double pv_fdr;

 public:

	/**
	 * \brief General constructor
	 * \params vector of field read from a tabulated file
	 */
	GenePvalue(const vector<string> & fields)
 	{
		gene = fields[0];
		diff = fields[1];
		caseSample = fields[2];
		controlSample = fields[3];
		fold_change = fields[4];
		hits = fields[5];
		nPermutation = fields[6];
		pv = atof(fields[7].c_str());
		pv_fdr = -1;
	}


	/**
	 * \brief return a string representation from all class members
	 */
	inline string toString()
	{
		return gene + "\t" + diff + "\t" + caseSample + "\t" + controlSample + "\t" + fold_change + "\t" + hits + "\t" + nPermutation + "\t" + numberToStr(pv) + "\t" + numberToStr(pv_fdr);
	}

	//> Get gene pvalue
	double getPv(){ return pv;}

	//> Get gene pvalue
	inline void setPvAdjusted(double newValue){ pv_fdr = newValue;}

};

/**
 * \brief splits a string a save its field in a vector
 * \param line to be splitted
 * \param vector <str> vecotr of string containing each line field
 */
inline void split(const string &line,vector <string> & strs)
{
	istringstream iss(line);
	copy(istream_iterator<string>(iss),istream_iterator<string>(),back_inserter(strs));
}

/**
 * \brief splits a string by tabs a save its field in a vector
 * \param line to be splitted
 * \param vector <str> vecotr of string containing each line field
 */
inline void splitTabs(const string &line,vector <string> & strs)
{
	istringstream iss(line);
	string token;
	while(std::getline(iss, token, '\t'))
	{
		strs.push_back(token);
	}
}



/**
 * \brief Parsers input files and creates gene structure
 * \param string nameFile Name of the file
 * \param vector<Gene> vector of gene information
 * \param string condition1 First condition group
 * \param string condition2 Secind condition group
 */
bool loadFileInfo(const string & nameFile,vector <Gene> & vGenes,const string & condition1,const string & condition2)
{
	ifstream geneExpresionFile;
	geneExpresionFile.open(nameFile.c_str());

	string lineRead = "";

	bool isHeaderRead = false;

	if (geneExpresionFile.is_open())
	{
		if (geneExpresionFile.good())
		{
			getline (geneExpresionFile,lineRead);
		}
		else
		{
			return false;
		}

		while (!lineRead.empty())
		{
			vector <string> vFields;
			split(lineRead,vFields);

			if (vFields.size() < 5)
			{
				break;
			}

			if(!isHeaderRead)
			{
				//READ FILE HEADER
				for (vector<string>::const_iterator iter = vFields.begin() + 4;iter != vFields.end(); ++iter)
				{
					Gene newGene;
					newGene.geneName = (*iter);
					vGenes.push_back(newGene);
				}
				isHeaderRead = true;
			}
			else
			{
				//READ NEW SAMPLE
				//group definition
				bool isCase = true;
				if (vFields[1].compare(condition2) == 0)
				{
					isCase = false;
				}

				vector<Gene>::iterator iterGenes = vGenes.begin();

				for (vector<string>::const_iterator iter = vFields.begin() + 4;iter != vFields.end(); ++iter)
				{
					(*iterGenes).addSample(vFields[0],isCase,atof((*iter).c_str()));
					++ iterGenes;
				}
			}

			getline (geneExpresionFile,lineRead);
		}
	}
	else
	{
		return false;
	}
	return true;
}

/* TODO:
 * def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):
    """
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1])
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    n = float(pvalues.shape[0])
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in xrange(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues
 *
 *
 */

/**
 * \brief Correct P value by Benjamini-Hochberg (FDR)
 * \param vPvalues  Vector of original pValues calculated
 * \param newPvalues Reference to the vector of new pValues to be corrected
 */
void correct_pvalues_fdr(const vector <double> & vPvalues,vector <double> & newPvalues)
{
	unsigned int nElements = vPvalues.size();
	vector <double> newValues;
	vector < pair <double,int> > values;

	unsigned int i = 0;
	for (vector<double>::const_iterator iter = vPvalues.begin();iter != vPvalues.end(); ++iter)
	{
		values.push_back(make_pair((*iter),i));
		newPvalues.push_back(-1);
		i++;
	}

	sort(values.begin(), values.end());
	reverse(values.begin(), values.end());

	int rank = 0;
	i = 0;
	for (vector< pair <double,int> >::const_iterator iter = values.begin();iter != values.end(); ++iter)
	{
		rank = nElements - i;
		newValues.push_back( (double) ((double)nElements/(double)rank) * (*iter).first );
		i++;
	}

	for (unsigned int a = 0; a < nElements; a++)
	{
		if (newValues[a] < newValues[a+1])
		{
			newValues[a+1] = newValues[a];
		}
	}


	//for (unsigned int a = 0; a < nElements; a++)
	unsigned int a = 0;
	for (vector< pair <double,int> >::const_iterator iter = values.begin();iter != values.end(); ++iter)
	{
		unsigned int index = (*iter).second;
		newPvalues[index] = newValues[a];
		a ++;
	}

}

/**
 * \brief Print application help
 */
void printHelp()
{
	cout << "diffExprpermutation: Find differentially expressed genes from a set of control and cases samples using a permutation strategy." << endl;
	cout << "diffExprpermutation -m MODE MODE PARAMETERS" << endl;

	cout << "Inputs:" << endl;
	cout << " -m MODE         SPLIT Splits original input file in a set of chunks. The goal is to perform a parallel running of the permutations process." << endl;
	cout << "                 PERM Run the permutations over a set of genes and samples to calculate a p value per gene." << endl;
	cout << "                 FDR p values adjustment per Benjamini-Hochber method, also known as False Discovery rate." << endl;
	cout << "To get more help run the application with your preferred mode and -h parameter:" << endl;
	cout << "                 diffExprPermutation -m SPLIT -h" << endl;
	cout << "                 diffExprPermutation -m PERM -h" << endl;
	cout << "                 diffExprPermutation -m FDR -h" << endl;

}

/**
 * \brief Get Header from files
 * \param nameFile Name of the file
 * \return First line in a file
 */
inline string getFileHeader(const string & fileName)
{
	ifstream inStreamFile;
	inStreamFile.open(fileName.c_str());
	string lineRead = "";

	if (inStreamFile.is_open())
	{
		if (inStreamFile.good())
		{
			getline (inStreamFile,lineRead);
		}
	}

	inStreamFile.close();
	return lineRead;
}

/**
 * \brief Adds a Gene Info to be printed in file of samples and expression values
 * \param firstGene True in case it is the first gene to be inserted in the vector of strings False in case all samples were already inserted just append expression values
 * \param vector chunkContents vector to be printed in file, first line is the header next line will be sample gene expression values
 * \param Gene reference to the information to be printed
 * \param string condition1 First group condition
 * \param string condition2 Second group condition
 */
void addGeneContents(bool firstGene,vector <string> & chunckContents, const Gene & gene, const string & condition1, const string & condition2)
{
	//Add Gene
	chunckContents[0].append(" " + gene.geneName); //Gene Name

	string group =  condition1;
	vector <bool>::const_iterator itGroup = gene.vGroup.begin();

	unsigned int samplePosition = 1;
	for (vector<SampleGene>::const_iterator itSamples = gene.vGeneValues.begin(); itSamples != gene.vGeneValues.end(); ++itSamples)
	{
		if(firstGene)
		{
			if((*itGroup)) group = condition1;
			else group = condition2;

			chunckContents.push_back((*itSamples).sampleName + " " + group +" 11111111 0.111111111111111 " + numberToStr((*itSamples).expressionValue));
			++itGroup;
		}
		else
		{
			chunckContents[samplePosition].append(" " + numberToStr((*itSamples).expressionValue));
			samplePosition ++;
		}
	}
}

/**
 * \brief Print a vector to string to a file
 * \param vector   String list of lines to be printed in file
 * \param unsigned int genesPrinted   Index to be printed as a part of the output file name
 * \param string outputPrefix File name output prefix
 */
void printChunck(vector <string> & chunckContents,unsigned int genesPrinted,const string & outputPrefix)
{
	string currentChunkFile;
	std::ofstream ofsChunkFile;

	//Chunk File Name
	currentChunkFile = outputPrefix + "_" + numberToStr(genesPrinted);
	//Open file descriptor
	ofsChunkFile.open(currentChunkFile.c_str(),std::ofstream::out);

	for (vector<string>::const_iterator it = chunckContents.begin(); it != chunckContents.end(); ++it)
	{
		ofsChunkFile << (*it) << endl;
	}

	ofsChunkFile.close();
	chunckContents.clear();
}

/**
 * \brief Load in a vector of GenePvalue objects a set of files located in folder
 * \param vector genePvalues Array of GenePvalue objects to be filled
 * \param string directory Folder path where are located the gene pvalues files
 * \param string header File Header
 */
void loadGenePValuesFromFolder(vector <GenePvalue> & genePvalues,const string & directory, string & header)
{
    DIR *dp;
	dp = opendir(directory.c_str());
	struct dirent *dirp;
	struct stat filestat;
	string filePath,lineRead;
	ifstream fin;

	if (dp == NULL)
	{
		cout << "Sorry!! Can not open directory:" << directory << endl;
		return;
	}

	while ((dirp = readdir( dp )))
	{
		filePath = directory + "/" + dirp->d_name;
		//3.1 If the file is a directory (or is in some way invalid) we'll skip it
		if (stat( filePath.c_str(), &filestat )) continue;
		if (S_ISDIR( filestat.st_mode ))         continue;
		//3.2 Read File
		fin.open( filePath.c_str() );

		bool isHeader = true;

		if (fin.is_open())
		{
			if (fin.good())
			{
				getline (fin,lineRead);
			}

			while (!lineRead.empty())
			{
				vector <string> vFields;
				splitTabs(lineRead,vFields);

				if(isHeader)
				{
					header = lineRead;
					isHeader = false;
				}
				else if(vFields.size() > 8)
				{
					genePvalues.push_back(GenePvalue(vFields));
				}

				getline (fin,lineRead);
			}
		}

		fin.close();
	}

	closedir( dp );
}

/**
 * \brief Print a vector of GenePvalue objects to a file
 * \param vector   String list of lines to be printed in file
 * \param unsigned int genesPrinted   Index to be printed as a part of the output file name
 * \param string outputPrefix File name output prefix
 */
void printCorrectedPvalues(vector <GenePvalue> & genesAdjusted,const string & header, const string & fileName)
{
	std::ofstream ofsFile;

	//Open file descriptor
	ofsFile.open(fileName.c_str(),std::ofstream::out);

	ofsFile << header << endl;

	for (vector<GenePvalue>::iterator it = genesAdjusted.begin(); it != genesAdjusted.end(); ++it)
	{
		ofsFile << (*it).toString() << endl;
	}

	ofsFile.close();
}

#endif /* COMMON_H_ */
