/**
 *  @file ReactionPath.cpp
 *  Implementation file for classes used in reaction path analysis.
 */

/*
 * $Author: dggoodwin $
 * $Revision: 1.1 $
 * $Date: 2007/05/04 14:27:23 $
 */

// Copyright 2001  California Institute of Technology

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ReactionPath.h"
#include "Kinetics.h"
#include "reaction_defs.h"
#include "Group.h"

using namespace std;

namespace Cantera {
     
    /// add a path to or from this node
    void SpeciesNode::addPath(Path* path) { 
        m_paths.push_back(path);
        if (path->begin() == this) m_out += path->flow();
        else if (path->end() == this) m_in += path->flow();
        else throw CanteraError("addPath","path added to wrong node");
    }
        
    void SpeciesNode::printPaths() {
        for (int i = 0; i < int(m_paths.size()); i++) {
            cout << m_paths[i]->begin()->name << " -->  "
                 << m_paths[i]->end()->name << ":   "
                 << m_paths[i]->flow() << endl;
        }
    }


    /** 
     * Construct a path connecting two species nodes.
     */   
    Path::Path(SpeciesNode* begin, SpeciesNode* end) 
        : m_a(begin), m_b(end), m_total(0.0) 
    { 
        begin->addPath(this);
        end->addPath(this);
    }


    /**
     * add a reaction to the path. Increment the flow from this
     * reaction, the total flow, and the flow associated with this
     * label.
     */
    void Path::addReaction(int rxnNumber, doublereal value, 
        string label) {
        m_rxn[rxnNumber] += value;
        m_total += value;
        if (label != "") m_label[label] += value;
    }


    /**
     * Write the label for a path connecting two species, indicating 
     * the percent of the total flow due to each reaction.
     */
    void Path::writeLabel(ostream& s, doublereal threshold) 
    {
        int nn = static_cast<int>(m_label.size());
        if (nn == 0) return;
        doublereal v;
        map<string, doublereal>::const_iterator i = m_label.begin();
        for (; i != m_label.end(); ++i) {
            v = i->second/m_total;
            if (nn == 1) s << i->first << "\\l";
            else if (v > threshold) {
                s << i->first;
                int percent = int(100*v + 0.5);
                if (percent < 100) 
                    s << " (" << percent << "%)\\l";
                else 
                    s << "\\l";
            }
        }
    }


    /**
     * Default constructor. 
     */
    ReactionPathDiagram::ReactionPathDiagram() {
        name = "reaction_paths";
        m_flxmax = 0.0;
        bold_color = "blue";
        normal_color = "steelblue";
        dashed_color = "gray";
        dot_options = "center=1;";
        m_font = RXNPATH_FONT;
        bold_min = 0.2;
        dashed_max = 0.0;
        label_min = 0.0;
        threshold = 0.005;
        flow_type = NetFlow;
        scale = -1;
        x_size = -1.0;
        y_size = -1.0;
        arrow_width = -5.0;
        show_details = false;
        arrow_hue = 0.6666;
        title = "";
        m_local = -1;
    }
        

    /**
     * Destructor. Deletes all nodes and paths in the diagram.
     */ 
    ReactionPathDiagram::~ReactionPathDiagram() 
    {
        // delete the nodes
        map<int, SpeciesNode*>::const_iterator i = m_nodes.begin();
        for (; i != m_nodes.end(); ++i) delete i->second;
        
        // delete the paths
        int nn = nPaths();
        int n;
        for (n = 0; n < nn; n++) delete m_pathlist[n];
    }


    vector_int ReactionPathDiagram::reactions() {
        int i, npaths = nPaths();
        double flmax = 0.0, flxratio;
        Path* p;
        for (i = 0; i < npaths; i++) 
        {
            p = path(i);
            if (p->flow() > flmax) flmax = p->flow();
        }
        m_rxns.clear();
        for (i = 0; i < npaths; i++) 
        {
            p = path(i);
            const Path::rxn_path_map& rxns = p->reactionMap();
            Path::rxn_path_map::const_iterator m = rxns.begin();
            for (; m != rxns.end(); ++m) {
                flxratio = m->second/flmax;
                if (flxratio > threshold) {
                    m_rxns[m->first] = 1;
                }
            }
        }
        vector_int r;
        map<int, int>::const_iterator begin = m_rxns.begin();
        for (; begin != m_rxns.end(); ++begin) r.push_back(abs(begin->first));
        return r;
    }

    void ReactionPathDiagram::add(ReactionPathDiagram& d) {
//        double f1, f2;
//         int nnodes = nNodes();
//         if (nnodes != d.nNodes()) {
//             throw CanteraError("ReactionPathDiagram::add",
//                 "number of nodes must be the same");
//         }
        int np = nPaths();
        int n, k1, k2;
        Path* p = 0;
        for (n = 0; n < np; n++) {
            p = path(n);
            k1 = p->begin()->number;
            k2 = p->end()->number;
            p->setFlow(p->flow() + d.flow(k1,k2));
        }
    }

    void ReactionPathDiagram::findMajorPaths(doublereal athreshold, int lda, 
        doublereal* a) {
        int nn = nNodes();
        int n, m, k1, k2;
        doublereal fl, netmax = 0.0;
        for (n = 0; n < nn; n++) {
            for (m = n+1; m < nn; m++) {
                k1 = m_speciesNumber[n];
                k2 = m_speciesNumber[m];
                fl = fabs(netFlow(k1,k2));
                if (fl > netmax) netmax = fl;
            }
        }
        for (n = 0; n < nn; n++) {
            for (m = n+1; m < nn; m++) {
                k1 = m_speciesNumber[n];
                k2 = m_speciesNumber[m];
                fl = fabs(netFlow(k1,k2));
                if (fl > athreshold*netmax) 
                    a[lda*k1 + k2] = 1;
            }
        }
    }

    void ReactionPathDiagram::writeData(ostream& s) {
        double f1, f2;
        int nnodes = nNodes();
        int i1, i2, k1, k2;
        s << title << endl;
        for (i1 = 0; i1 < nnodes; i1++) 
        {
            k1 = m_speciesNumber[i1];
            s << m_nodes[k1]->name << " ";
        }
        s << endl;
        for (i1 = 0; i1 < nnodes; i1++) 
        {
            k1 = m_speciesNumber[i1];
            for (i2 = i1+1; i2 < nnodes; i2++) 
            {
                k2 = m_speciesNumber[i2];
                f1 = flow(k1, k2);
                f2 = flow(k2, k1);
                //if (f1 > 0.001 || f2 > 0.001) {
                s << m_nodes[k1]->name << " " << m_nodes[k2]->name
                  << " " << f1 << " " << -f2 << endl;
                //}
            }
        }
    }


    /** 
     *  Export the reaction path diagram. This method writes to stream
     *  \c s the commands for the 'dot' program in the \c GraphViz
     *  package from AT&T. (GraphViz may be downloaded from
     *  www.graphviz.org.)
     *
     *  To generate a postscript reaction path diagram from the
     *  output of this method saved in file paths.dot, for example, give
     *  the command: 
     *  \code 
     *  dot -Tps paths.dot > paths.ps 
     *  \endcode
     *  To generate a GIF image, replace -Tps with -Tgif
     */    
    void ReactionPathDiagram::exportToDot(ostream& s) 
    {
        int i;
        doublereal flxratio, flmax = 0.0, lwidth;
        //s.flags(std::ios_base::showpoint+std::ios_base::fixed);
        s.precision(3);

        // a directed graph
        s << "digraph " << name << " {" << endl;

        // the graph will be no larger than x_size, y_size
        if (x_size > 0.0) 
        {
            if (y_size < 0.0) y_size = x_size;
            s << "size = \"" 
              << x_size << "," 
              << y_size << "\";" 
              << endl;
        } 

        //s << "color = white;" << endl;
        if (dot_options != "") 
            s << dot_options << endl;

        int npaths = nPaths();
        Path* p;

        int nnodes = nNodes();
        int kbegin, kend, i1, i2, k1, k2;
        double flx;

        // draw paths representing net flows
        if (flow_type == NetFlow) 
        {

            // if no scale was specified, normalize
            // net flows by the maximum net flow
            if (scale <= 0.0) 
            {
                for (i1 = 0; i1 < nnodes; i1++) 
                {
                    k1 = m_speciesNumber[i1];
                    node(k1)->visible  = false;
                    for (i2 = i1+1; i2 < nnodes; i2++) 
                    {
                        k2 = m_speciesNumber[i2];
                        flx = netFlow(k1, k2);
                        if (flx < 0.0) flx = -flx;
                        if (flx > flmax) flmax = flx;
                    }
                }
            }
            else 
                flmax = scale;
            
            if (flmax < 1.e-10) flmax = 1.e-10;

            // loop over all unique pairs of nodes

            for (i1 = 0; i1 < nnodes; i1++) 
            {
                k1 = m_speciesNumber[i1];
                for (i2 = i1+1; i2 < nnodes; i2++) 
                {
                    k2 = m_speciesNumber[i2];
                    flx = netFlow(k1, k2);
                    if (m_local >= 0) {
                        if (k1 != m_local && k2 != m_local) flx = 0.0;
                    }
                    if (flx != 0.0) 
                    {
                        // set beginning and end of the path based on the
                        // sign of the net flow

                        if (flx > 0.0) 
                        {
                            kbegin = k1;
                            kend = k2;
                            flxratio = flx/flmax;
                        }
                        else 
                        {
                            kbegin = k2;
                            kend = k1;
                            flxratio = -flx/flmax;
                        }

                        // write out path specification if the net flow
                        // is greater than the threshold

                        if (flxratio >= threshold) 
                        {
                            // make nodes visible
                            node(kbegin)->visible  = true;
                            node(kend)->visible    = true;

                            s << "s" << kbegin << " -> s" << kend;

                            if (arrow_width < 0) {
                                lwidth = 1.0 - 4.0 
                                         * log10(flxratio/threshold)/log10(threshold) + 1.0;
                                s <<  "[fontname=\""+m_font+"\", style=\"setlinewidth("
                                  << lwidth << ")\"";
                                s << ", arrowsize=" 
                                  <<  min(6.0, 0.5*lwidth); 
                            } 
                            else 
                            {
                                s <<  ", style=\"setlinewidth(" 
                                  <<  arrow_width << ")\"";
                                s << ", arrowsize=" << flxratio + 1; 
                            }

                            doublereal hue = 0.7;
                            doublereal bright = 0.9;
                            s << ", color=" << "\"" << hue << ", " 
                              << flxratio + 0.5
                              << ", " << bright << "\"" << endl;

                            if (flxratio > label_min) 
                            {
                                s << ", label=\" " << flxratio;
                                if (show_details) {
                                    if (flow(kbegin, kend) > 0.0) {
                                        s << "\\l fwd: " 
                                          << flow(kbegin, kend)/flmax << "\\l";
                                        path(kbegin, kend)->writeLabel(s);
                                    }
                                    if (flow(kend, kbegin) > 0.0) {
                                        s << " \\l rev: " 
                                          << flow(kend,kbegin)/flmax << "\\l";
                                        path(kend, kbegin)->writeLabel(s);
                                    }
                                }
                                s << "\"";
                            }
                            s << "];" << endl;
                        }
                    }
                }
            }
        }

        else {
            for (i = 0; i < npaths; i++) 
            {
                p = path(i);
                if (p->flow() > flmax) flmax = p->flow();
            }

            for (i = 0; i < npaths; i++) {
                p = path(i);
                flxratio = p->flow()/flmax;
                if (m_local >= 0) {
                    if (p->begin()->number != m_local 
                        && p->end()->number != m_local) flxratio = 0.0;
                }
                if (flxratio > threshold) {
                    p->begin()->visible = true;       
                    p->end()->visible = true;
                    s << "s" << p->begin()->number 
                        << " -> s" << p->end()->number; 

                    if (arrow_width < 0) {
                        lwidth = 1.0 - 4.0 * log10(flxratio/threshold)/log10(threshold)
 + 1.0;
                        s <<  "[fontname=\""+m_font+"\", style=\"setlinewidth("
                            //<< 1.0 - arrow_width*flxratio
                          << lwidth
                          << ")\"";
                        s << ", arrowsize=" 
                          <<  min(6.0, 0.5*lwidth); // 1 - arrow_width*flxratio; 
                    }
                    else 
                    {
                        s <<  ", style=\"setlinewidth(" 
                          <<  arrow_width << ")\"";
                        s << ", arrowsize=" << flxratio + 1; 
                    }
                    doublereal hue = 0.7; //2.0/(1.0 + pow(log10(flxratio),2)) ;
                    doublereal bright = 0.9;
                    s << ", color=" << "\"" << hue << ", " << flxratio + 0.5
                      << ", " << bright << "\"" << endl;

                    if (flxratio > label_min) {
                        s << ", label = \" " << flxratio;
                        if (show_details) {
                            s << "\\l"; p->writeLabel(s);
                        }
                        s << "\"";
                    }
                    s << "];" << endl;
                }
            }
        }
        s.precision(2);
        map<int, SpeciesNode*>::const_iterator b = m_nodes.begin();
        for (; b != m_nodes.end(); ++b) {
            if (b->second->visible) {
            s << "s" << b->first << " [ fontname=\""+m_font+"\", label=\"" << b->second->name
              //<< " \\n " << b->second->value
              << "\"];" << endl;
            }
        }
        s << " label = " << "\"" << "Scale = " 
          << flmax << "\\l " << title << "\";" << endl; //created with Cantera (www.cantera.org)\\l\";"
        s  << " fontname = \""+m_font+"\";" << endl << "}" << endl;
    }


    void ReactionPathDiagram::addNode(int k, string nm, doublereal x) {
        if (!m_nodes[k]) {
            m_nodes[k] = new SpeciesNode;
            m_nodes[k]->number = k;
            m_nodes[k]->name = nm;
            m_nodes[k]->value = x;
            m_speciesNumber.push_back(k);
        }
    }

    void ReactionPathDiagram::linkNodes(int k1, int k2, int rxn, 
        doublereal value, string legend) {
        SpeciesNode* begin = m_nodes[k1];
        SpeciesNode* end   = m_nodes[k2];
        Path* ff = m_paths[k1][k2];
        if (!ff) {
            ff= new Path(begin, end);
            m_paths[k1][k2] = ff;
            m_pathlist.push_back(ff);
        }
        ff->addReaction(rxn, value, legend);
        m_rxns[rxn] = 1;
        if (ff->flow() > m_flxmax) m_flxmax = ff->flow();
    }
    
    vector_int ReactionPathDiagram::species(){
        return m_speciesNumber;
    }


    /**
     *  analyze a reaction to determine which reactants lead to which products.
     */
    int ReactionPathBuilder::findGroups(ostream& logfile, Kinetics& s) 
    {
        m_groups.resize(m_nr);
	map<int, int> net;

        for (int i = 0; i < m_nr; i++)             // loop over reactions 
        {
            logfile << endl << "Reaction " << i+1 << ": " 
                    << s.reactionString(i);

	    int nrnet = m_reac[i].size();
	    int npnet = m_prod[i].size();
            const vector_int& r = s.reactants(i);
            const vector_int& p = s.products(i);

            int nr = s.reactants(i).size();
            int np = s.products(i).size();
 
            Group b0, b1, bb;

            vector<string>& e = m_elementSymbols;

            const vector<grouplist_t>& rgroups = s.reactantGroups(i);
            const vector<grouplist_t>& pgroups = s.productGroups(i);
            

            if (m_determinate[i]) {
                logfile << " ... OK." << endl;
            }

            else if (rgroups.size() > 0) {
                logfile << " ... specified groups." << endl;
                int nrg = static_cast<int>(rgroups.size());
                int npg = static_cast<int>(pgroups.size());
                int kr, kp, ngrpr, ngrpp;
                Group gr, gp;

                if (nrg != nr || npg != np) return -1; 
                
                // loop over reactants 
                for (int igr = 0; igr < nrg; igr++) {
                    kr = r[igr];
                    ngrpr = static_cast<int>(rgroups[igr].size());

                    // loop over products
                    for (int igp = 0; igp < npg; igp++) {
                        kp = p[igp];
                        ngrpp = static_cast<int>(pgroups[igp].size());

                        // loop over pairs of reactant and product groups
                        for (int kgr = 0; kgr < ngrpr; kgr++) {
                            gr = Group(rgroups[igr][kgr]);
                            for (int kgp = 0; kgp < ngrpp; kgp++) {
                                gp = Group(pgroups[igp][kgp]);
                                if (gr == gp) {
                                    m_transfer[i][kr][kp] = gr;
                                }
                            }
                        }
                    }
                }
            }

            else if (nrnet == 2 && npnet == 2) 
            {
                // indices for the two reactants
                int kr0 = m_reac[i][0];
                int kr1 = m_reac[i][1];

                // indices for the two products
                int kp0 = m_prod[i][0];
                int kp1 = m_prod[i][1];

                // references to the Group objects representing the
                // reactants 
                const Group& r0 = m_sgroup[kr0];
                const Group& r1 = m_sgroup[kr1];
                const Group& p0 = m_sgroup[kp0];
                const Group& p1 = m_sgroup[kp1];

                const Group *group_a0=0, *group_b0=0, *group_c0=0, 
                            *group_a1=0, *group_b1=0, *group_c1=0;
                b0 = p0 - r0;
                b1 = p1 - r0;
                if (b0.valid() && b1.valid()) {
                    logfile << " ... ambiguous." << endl;
                }
                else if (!b0.valid() && !b1.valid()) {
                    logfile << " ... cannot express as A + BC = AB + C" << endl;
                }
                else logfile << endl;

                if (b0.valid()) {
                    if (b0.sign() > 0) {
                        group_a0 = &r0;
                        group_b0 = &b0;
                        group_c0 = &p1;
                        m_transfer[i][kr0][kp0] = r0;
                        m_transfer[i][kr1][kp0] = b0;
                        m_transfer[i][kr1][kp1] = p1;
                    }
                    else {
                        group_a0 = &r1;
                        group_c0 = &p0;
                        b0 *= -1;
                        group_b0 = &b0;
                        m_transfer[i][kr1][kp1] = r1;
                        m_transfer[i][kr0][kp1] = b0;
                        m_transfer[i][kr0][kp0] = p0;
                    }
                    logfile << "     "; 
                    group_a0->fmt(logfile,e); 
                    logfile << " + "; 
                    group_b0->fmt(logfile,e); 
                    group_c0->fmt(logfile,e); 
                    logfile << " = ";
                    group_a0->fmt(logfile,e); 
                    group_b0->fmt(logfile,e); 
                    logfile << " + ";  
                    group_c0->fmt(logfile,e); 
                    if (b1.valid()) 
                        logfile << "   [<= default] " << endl;
                    else 
                        logfile << endl;
                }


                if (b1.valid()) {
                    if (b1.sign() > 0) {
                        group_a1 = &r0;
                        group_b1 = &b1;
                        group_c1 = &p0;
                        if (!b0.valid()) {
                            m_transfer[i][kr0][kp1] = r0;
                            m_transfer[i][kr1][kp1] = b0;
                            m_transfer[i][kr1][kp0] = p0;
                        }
                    }
                    else {
                        group_a1 = &r1;
                        group_c1 = &p1;
                        b1 *= -1;
                        group_b1 = &b1;
                        if (!b0.valid()) {
                            m_transfer[i][kr1][kp0] = r1;
                            m_transfer[i][kr0][kp0] = b0;
                            m_transfer[i][kr0][kp1] = p1;
                        }
                    }
                    logfile << "     ";
                    group_a1->fmt(logfile,e); logfile << " + "; 
                    group_b1->fmt(logfile,e); 
                    group_c1->fmt(logfile,e); logfile << " = ";
                    group_a1->fmt(logfile,e); group_b1->fmt(logfile,e); 
                    logfile << " + ";  group_c1->fmt(logfile,e); 
                    logfile << endl;
                }
            }
            else {
                logfile << "... cannot parse. [ignored]" << endl;
            }
        }
        return 1;
    }
                        
    void ReactionPathBuilder::writeGroup(ostream& out, const Group& g) 
    {
        g.fmt(out, m_elementSymbols);
    }

    void ReactionPathBuilder::findElements(Kinetics& kin) {

        string ename;
        m_enamemap.clear();
        m_nel = 0;
        int i, np = kin.nPhases();
        ThermoPhase* p;
        map<string, int> enamemap;
        for (i = 0; i < np; i++) {
            p = &kin.thermo(i);
            // iterate over the elements in this phase
            int m, nel = p->nElements();
            for (m = 0; m < nel; m++) {
                ename = p->elementName(m);

                // if no entry is found for this element name, then
                // it is a new element. In this case, add the name
                // to the list of names, increment the element count, 
                // and add an entry to the name->(index+1) map.
                if (m_enamemap.find(ename) == m_enamemap.end()) {
                    m_enamemap[ename] = m_nel + 1;
                    m_elementSymbols.push_back(ename);
                    m_nel++;
                }
            }
        }
        m_atoms.resize(kin.nTotalSpecies(), m_nel, 0.0);
        string sym;
        int k, ip, nsp, mlocal, kp, m;
        // iterate over the elements
        for (m = 0; m < m_nel; m++) {
            sym = m_elementSymbols[m];
            k = 0;
            // iterate over the phases
            for (ip = 0; ip < np; ip++) {
                phase_t* p = &kin.thermo(ip);
                nsp = p->nSpecies();
                mlocal = p->elementIndex(sym);    
                for (kp = 0; kp < nsp; kp++) {
                    if (mlocal >= 0) {
                        m_atoms(k, m) = p->nAtoms(kp, mlocal);
                    }
                    k++;
                }
            }
        }
    }



    int ReactionPathBuilder::init(ostream& logfile, Kinetics& kin) {
        //m_warn.clear();
        m_transfer.clear();

        //const Kinetics::thermo_t& ph = kin.thermo();

        m_elementSymbols.clear();
        findElements(kin);
        //m_nel = ph.nElements();
        m_ns = kin.nTotalSpecies(); //ph.nSpecies();
        m_nr = kin.nReactions();

        int m, i;
        //for (m = 0; m < m_nel; m++) {
        //    m_elementSymbols.push_back(ph.elementName(m));
        //}

        // all reactants / products, even ones appearing on both sides
        // of the reaction
        // mod 8/18/01 dgg
        vector<vector_int > allProducts;
        vector<vector_int > allReactants;
        for (i = 0; i < m_nr; i++) {
            allReactants.push_back(kin.reactants(i));
            allProducts.push_back(kin.products(i));
        }

        // m_reac and m_prod exclude indices for species that appear on
        // both sides of the reaction, so that the diagram contains no loops.

        m_reac.resize(m_nr);
        m_prod.resize(m_nr);

        m_ropf.resize(m_nr);
        m_ropr.resize(m_nr);
        m_determinate.resize(m_nr);

        m_x.resize(m_ns);  // not currently used ?
        m_elatoms.resize(m_nel, m_nr);
        
        int nr, np, n, k;
        int nmol;
        map<int, int> net;

        for (i = 0; i < m_nr; i++) {

            // construct the lists of reactant and product indices, not 
            // including molecules that appear on both sides. 

	    m_reac[i].clear();
	    m_prod[i].clear();
	    net.clear();
            nr = allReactants[i].size();
            np = allProducts[i].size();
	    for (int ir = 0; ir < nr; ir++) net[allReactants[i][ir]]--;
	    for (int ip = 0; ip < np; ip++) net[allProducts[i][ip]]++;

	    for (k = 0; k < m_ns; k++) {
	      if (net[k] < 0) {
                  nmol = -net[k];
                  for (int jr = 0; jr < nmol; jr++) m_reac[i].push_back(k);
	      }
	      else if (net[k] > 0) {
                  nmol = net[k];
                  for (int jp = 0; jp < nmol; jp++) m_prod[i].push_back(k);
	      }
	    }

 	    int nrnet = m_reac[i].size();
            // 	    int npnet = m_prod[i].size();

            // compute number of atoms of each element in each reaction,
            // excluding molecules that appear on both sides of the
            // reaction. We only need to compute this for the reactants,
            // since the elements are conserved.

            for (n = 0; n < nrnet; n++) {
                k = m_reac[i][n];
                for (int m = 0; m < m_nel; m++) {
                    m_elatoms(m,i) += m_atoms(k,m); //ph.nAtoms(k,m);
                }
            }            
        }

        // build species groups
        vector_int comp(m_nel);
        m_sgroup.resize(m_ns);
        int j;
        for (j = 0; j < m_ns; j++) {
            for (int m = 0; m < m_nel; m++) comp[m] = int(m_atoms(j,m)); //ph.nAtoms(j,m));
            m_sgroup[j] = Group(comp);
        }


        // determine whether or not the reaction is "determinate", meaning
        // that there is no ambiguity about which reactant is the source for 
        // any element in any product. This is false if more than one
        // reactant contains a given element, *and* more than one product
        // contains the element. In this case, additional information is 
        // needed to determine the partitioning of the reactant atoms of 
        // that element among the products.
 
        int nar, nap;
        for (i = 0; i < m_nr; i++) {
            nr = m_reac[i].size();
            np = m_prod[i].size();
            m_determinate[i] = true;
            for (m = 0; m < m_nel; m++) {
                nar = 0;
                nap = 0;
                for (j = 0;  j < nr; j++) {
                    //                    if (ph.nAtoms(m_reac[i][j],m) > 0) nar++;
                    if (m_atoms(m_reac[i][j],m) > 0) nar++;
                }
                for (j = 0;  j < np; j++) {
                    if (m_atoms(m_prod[i][j],m) > 0) nap++;
                }
                if (nar > 1 && nap > 1) {
                    m_determinate[i] = false; break;
                }
            }
        }

        findGroups(logfile, kin);
        return 1;
    }

    string reactionLabel(int i, int kr, int nr, const vector_int& slist, 
        const Kinetics& s) {
        
        //int np = s.nPhases();
        string label = "";
        int l;
        for (l = 0; l < nr; l++) {
            if (l != kr) 
                label += " + "+ s.kineticsSpeciesName(slist[l]);
        }
        if (s.reactionType(i) == THREE_BODY_RXN)
            label += " + M "; 
        else if (s.reactionType(i) == FALLOFF_RXN)
            label += " (+ M)";
        return label;
    }


    int ReactionPathBuilder::build(Kinetics& s, 
        string element, ostream& output, ReactionPathDiagram& r, bool quiet) 
    {    
        int i, nr, np, kr, kp, kkr, kkp;
        doublereal f, ropf, ropr, fwd, rev;
        string fwdlabel, revlabel;
        map<int, int> warn;

        doublereal threshold = 0.0;
        bool fwd_incl, rev_incl, force_incl;

        //        const Kinetics::thermo_t& ph = s.thermo();
        int m = m_enamemap[element]-1; //ph.elementIndex(element);

        r.element = element;
        if (m < 0) return -1;
        
        //int k;
        int kk = s.nTotalSpecies();

        s.getFwdRatesOfProgress(DATA_PTR(m_ropf));
        s.getRevRatesOfProgress(DATA_PTR(m_ropr));

        //ph.getMoleFractions(m_x.begin());

        //doublereal sum = 0.0;
        //for (k = 0; k < kk; k++) {
        //    sum += m_x[k] * ph.nAtoms(k,m);
        //}
        //sum *= ph.molarDensity();

        // species explicitly included or excluded
        vector<string>& in_nodes = r.included();
        vector<string>& out_nodes = r.excluded();
        int nin = static_cast<int>(in_nodes.size());
        int nout = static_cast<int>(out_nodes.size());

        vector_int status;
	status.resize(kk,0);
        for (int ni = 0; ni < nin; ni++) 
            status[s.kineticsSpeciesIndex(in_nodes[ni])] = 1;
        for (int ne = 0; ne < nout; ne++) 
            status[s.kineticsSpeciesIndex(out_nodes[ne])] = -1;

        for (i = 0; i < m_nr; i++) 
        {
            ropf = m_ropf[i];
            ropr = m_ropr[i];

            // loop over reactions involving element m
            if (m_elatoms(m, i) > 0) 
            {
                nr = m_reac[i].size();
                np = m_prod[i].size();

                for (kr = 0; kr < nr; kr++) 
                {
                    kkr = m_reac[i][kr];
                    int l;

                    fwdlabel = reactionLabel(i, kr, nr, m_reac[i], s);
 
                    for (kp = 0; kp < np; kp++) 
                    {
                        kkp = m_prod[i][kp];
                        revlabel = "";
                        for (l = 0; l < np; l++) {
                            if (l != kp) 
                                revlabel += " + "+ s.kineticsSpeciesName(m_prod[i][l]);
                        }
                        if (s.reactionType(i) == THREE_BODY_RXN)
                            revlabel += " + M "; 
                        else if (s.reactionType(i) == FALLOFF_RXN)
                            revlabel += " (+ M)"; 


                        // calculate the flow only for pairs that are
                        // not the same species, both contain atoms of 
                        // element m, and both are allowed to appear in
                        // the diagram

                        if ((kkr != kkp) && (m_atoms(kkr,m) > 0 
                            &&  m_atoms(kkp,m) > 0)
                            && status[kkr] >= 0 && status[kkp] >= 0) 
                        {

                            // if neither species contains the full
                            // number of atoms of element m in the
                            // reaction, then we must consider the
                            // type of reaction to determine which
                            // reactant species was the source of a
                            // given m-atom in the product
 
                            if ( (m_atoms(kkp,m) < m_elatoms(m, i)) && 
                                (m_atoms(kkr,m) < m_elatoms(m, i)) ) 
                            {
                                map<int, map<int, Group> >& g = m_transfer[i];
                                if (g.empty()) {
                                    if (!warn[i]) {
                                        if (!quiet) {
                                            output << endl;
                                            output << "*************** REACTION IGNORED ***************" << endl; 
                                            output << "Warning: no rule to determine partitioning of " << element
                                                   << endl << " in reaction " << s.reactionString(i) << "." << endl 
                                                   << "*************** REACTION IGNORED **************" << endl;
                                            output << endl;
                                            warn[i] = 1;
                                        }
                                    }                                    
                                    f = 0.0;
                                }
                                else {
                                    if (!g[kkr][kkp]) f = 0.0;
                                    else f = g[kkr][kkp].nAtoms(m);
                                }
                            }

                            // no ambiguity about where the m-atoms come
                            // from or go to. Either all reactant m atoms
                            // end up in one product, or only one reactant 
                            // contains all the m-atoms. In either case,
                            // the number of atoms transferred is given by
                            // the same expression.

                            else {
                                f = m_atoms(kkp,m) * m_atoms(kkr,m) / m_elatoms(m, i);
                            }

                            fwd = ropf*f;
                            rev = ropr*f;
                            force_incl = ((status[kkr] == 1) || (status[kkp] == 1));

                            fwd_incl = ((fwd > threshold) || 
                                (fwd > 0.0 && force_incl));
                            rev_incl = ((rev > threshold) || 
                                (rev > 0.0 && force_incl));
                            if (fwd_incl || rev_incl) 
                            {
                                if (!r.hasNode(kkr)) {
                                    r.addNode(kkr, s.kineticsSpeciesName(kkr), m_x[kkr]);
                                }
                                if (!r.hasNode(kkp)) {
                                    r.addNode(kkp, s.kineticsSpeciesName(kkp), m_x[kkp]);
                                }
                            }
                            if (fwd_incl) {
                                r.linkNodes(kkr, kkp, i, fwd, fwdlabel);
                            }
                            if (rev_incl) {
                                r.linkNodes(kkp, kkr, -i, rev, revlabel);
                            }		
                        }
                    }
                }
            }
        }
        return 1;
    }


}
