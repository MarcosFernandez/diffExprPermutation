/*
 * main.cpp
 *
 *  Created on: 30 Xan, 2015
 *      Author: marcos
 */
#include "common.h"

extern "C" {
#include "perm.h"
}


/**
 * \brief Print application help
 */
void printHelp() {
  cout << "diffExprpermutation: Find differentially expressed genes from a set of control and cases samples using a permutation strategy." << endl;
  cout << "diffExprPermutation -f input [--c1 condition1 --c2 condition2 -n permutations -H stopHits -s statistic] -o outFile" << endl;
  cout << "Inputs:" << endl;
  cout << " -f input        Space separated table. Format: sampleName group lib.size norm.factors gene1 gene2 ... geneN" << endl;
  cout << "Outputs:" << endl;
  cout << " -o outFile      Output file name" << endl;
  cout << "Options:" << endl;
  cout << " --c1	condition1            Condition that determine one of two groups [default: case]" << endl;
  cout << " --c2	condition2            Condition that determine other group [default: control]" << endl;
  cout << " -s statistic                Statistic to compute pvalue median|perc25|perc75|x [default: median]" << endl;
  cout << " -p percentile mode          Mode for selection of percentile auto|linear|nearest [default: auto]" << endl;
  cout << " -t n_threads                Number of threads [default: 1]" << endl;
}


/**
 * \brief Check Arguments
 * \param string fileInput - Name input file
 * \param string outFile -   Name output file
 * \param unsigned int chunks - Number of chunks to create
 * \param string condition1 - First condition group. Usually case.
 * \param string condition1 - Second condition group. Usually control.
 */
inline bool argumentChecking(const string &fileInput,
                                     const string &outFile,
                                     const string &condition1,
                                     const string &condition2) 
{
  bool bWrong = false;

  if (fileInput.empty()) 
  {
    cout << "Sorry!! No input file was specified!!" << endl;
    return true;
  }

  if (outFile.empty()) 
  {
    cout << "Sorry!! No output file was specified!!" << endl;
    return true;
  }

  if (condition1.empty()) 
  {
    cout << "Sorry!! Condition group 1 is empty!!" << endl;
    return true;
  }

  if (condition2.empty()) 
  {
    cout << "Sorry!! Condition group 2 is empty!!" << endl;
    return true;
  }

  return bWrong;
}


int main(int argc, char **argv) 
{
  string fileInput = "";
  string outFile = "";
  string condition1 = "case";
  string condition2 = "control";
  string percentile_mode = "auto";
  cp_mode pc_mode = AUTO;
  int n_threads = 1;
  string statistic = "median";
  double fStatisticValue = 0;
  bool doMedian = true;
  vector<Gene> vGenes; // vector of genes where each gene has a vector of sampleGenes, each sampleGene contains sample name expression value and group
  /**
   * BRACA1 -> A,true,0.75
   *        -> B,false,0.85
   *        ...
   * BRACA2 -> A,true,0.15
   *        -> B,false,0.20
   *        ...
   */

  // 1.Process parameters
  for (int i = 1; i < argc; i++) 
  {
    if (strcmp(argv[i], "-f") == 0)
    {
        fileInput = argv[++i];
    } 
    else if (strcmp(argv[i], "-o") == 0) 
    {
        outFile = argv[++i];
    } 
    else if (strcmp(argv[i], "--c1") == 0)
    {
        condition1 = argv[++i];
    }
    else if (strcmp(argv[i], "--c2") == 0)
    {
        condition2 = argv[++i];
    }
    else if (strcmp(argv[i], "-s") == 0)
    {
        statistic = argv[++i];
    }
    else if (strcmp(argv[i], "-p") == 0) 
    {
        percentile_mode = argv[++i];
    } 
    else if (strcmp(argv[i], "-t") == 0) 
    {
        n_threads = atoi(argv[++i]);
        if (n_threads < 1) n_threads = 1;
    }
    else if (strcmp(argv[i],"-h") == 0)
    {
	printHelp();
	return 0;
    }
  }

  // Check Arguments
  if(argumentChecking(fileInput, outFile, condition1, condition2)) 
  {
    return -1;
  }

  // Updates statistic
  string headerOutput = "gene\tdiff_median\tmedianCase\tmedianControl\tfold_change\tmedian_pv\tmedian_pv_fdr";

  if (statistic.compare("perc25") == 0) 
  {
    fStatisticValue = 25.0;
    doMedian = false;
    headerOutput = "gene\tdiff_lowerq\tlowerqCase\tlowerqControl\tfold_change\tlowerq_pv\tlowerq_pv_fdr";
  } 
  else if (statistic.compare("perc75") == 0) 
  {
    fStatisticValue = 75.0;
    doMedian = false;
    headerOutput = "gene\tdiff_UpQ\tupperqCase\tupperqControl\tfold_change\tupperq_pv\tupper_pv_fdr";
  } 
  else 
  {
    char *p;
    double x = strtod(statistic.c_str(), &p);
    if (x > 0.0 && x < 100.0) 
    {
      fStatisticValue = x;
      doMedian = false;
      ostringstream s;
      s << "gene\tdiff_" << x << "%\t" << x << "\%_Case\t" << x << "\%_Control\tfold_change\t" << x << "\%_pv\t" << x << "\%_pv_fdr";
      headerOutput = s.str();
    }
   }
  
   if (percentile_mode.compare("auto") == 0) pc_mode = AUTO;
   else if (percentile_mode.compare("linear") == 0) pc_mode = LINEAR_INTERPOLATION;
   else if (percentile_mode.compare("nearest") == 0) pc_mode = NEAREST_RANK;
   else 
   {
         cerr << "Percentile mode '" << percentile_mode << "' not recognized" << endl;
         return -1;
   }

   // Parsing Input file
   if (!loadFileInfo(fileInput, vGenes, condition1, condition2)) 
   {
       cerr << "Sorry!! Can not open file " << fileInput << endl;
       return -1;
   }

   // Allocate and make C structure for permutation routine
   struct perm_data pdata;
   pdata.n_cases = vGenes.begin()->nCases;
   pdata.n_controls = vGenes.begin()->nControls;
   int n_samples = pdata.n_cases + pdata.n_controls;
   pdata.n_var = vGenes.size();
   pdata.data = (float *)malloc(sizeof(float) * pdata.n_var * (4 + n_samples));
   pdata.fold_change = pdata.data + pdata.n_var * n_samples;
   pdata.pc_cases = pdata.fold_change + pdata.n_var;
   pdata.pc_controls = pdata.pc_cases + pdata.n_var;
   pdata.pc_diff = pdata.pc_controls + pdata.n_var;
   pdata.grp = (group *)malloc(sizeof(group) * n_samples);
   pdata.p_values = (double *)malloc(sizeof(double) * pdata.n_var);

   // Copy expression data
   float *tp = pdata.data;
   for (vector<Gene>::iterator gene_it = vGenes.begin(); gene_it != vGenes.end(); ++gene_it) 
   {
     for (vector<SampleGene>::const_iterator iter = gene_it->vGeneValues.begin();iter != gene_it->vGeneValues.end(); ++iter) 
     {
       (*tp++) = (float)iter->expressionValue;
     }
   }

   // Copy group information
   int ix = 0;
   for (vector<bool>::const_iterator iter = vGenes.begin()->vGroup.begin(); iter != vGenes.begin()->vGroup.end(); ++iter)
   {
       pdata.grp[ix++] = *iter ? CASE : CONTROL;
   }
  
   // Calculate exact permutation p-values
   check_percentile(doMedian ? 50.0 : fStatisticValue, pc_mode, n_threads,&pdata);

   vector<double> vPvalues;

   int i = 0;
   for (vector<Gene>::iterator iter = vGenes.begin(); iter != vGenes.end(); ++iter) 
   {
       // Copy results from pdata struture
       (*iter).originalDiff = (double)pdata.pc_diff[i];
       (*iter).originalMedianCases = (double)pdata.pc_cases[i];
       (*iter).originalMedianControl = (double)pdata.pc_controls[i];
       (*iter).foldChange = pdata.fold_change[i];
       (*iter).pValue = pdata.p_values[i];

       // Add pvalue to a vector of pvalues for being correcte by FDR
       vPvalues.push_back((*iter).pValue);
       i++;
   }

   vector<double> correctedPvalues;
   correct_pvalues_fdr(vPvalues, correctedPvalues);

   // Print to file
   std::ofstream outfile(outFile.c_str(), std::ofstream::out);

   // Header File
   outfile << headerOutput << endl;

   vector<double>::const_iterator iterCorrected = correctedPvalues.begin();
   outfile.precision(15);

   for (vector<Gene>::const_iterator iter = vGenes.begin(); iter != vGenes.end();++iter) 
   {
       outfile << (*iter).geneName << "\t" << (*iter).originalDiff << "\t" << (*iter).originalMedianCases;
       outfile << "\t" << (*iter).originalMedianControl << "\t" << (*iter).foldChange << "\t" << (*iter).pValue << "\t" << (*iterCorrected) << endl;
       ++iterCorrected;
   }

   free(pdata.grp);
   free(pdata.data);
   free(pdata.p_values);

   outfile.close();
 
  return 0;


}
