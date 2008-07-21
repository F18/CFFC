
#ifndef _FILTER_STATE_INCLUDED
#define _FILTER_STATE_INCLUDED

class Filter_State {
public:
    Filter_State(void){
        member = ZERO;
    }
    double member;
    int NumVar() {
        return 1;
    }
    
    double &operator[](int index);		
    const double &operator[](int index) const;
    friend ostream& operator << (ostream &out_file, const Filter_State &F);
    friend istream& operator >> (istream &in_file,  Filter_State &F); 
};

inline double& Filter_State::operator[](int index) {
    return member;
}

inline const double& Filter_State::operator[](int index) const {
    return member;
}


inline ostream& operator << (ostream &out_file, 
                             const Filter_State &F) {
    out_file.setf(ios::scientific);
    out_file << " " << F.member;
    out_file.unsetf(ios::scientific);
    return (out_file);
}
inline istream& operator >> (istream &in_file, 
                             Filter_State &F) {
    in_file.setf(ios::skipws);
    in_file >> F.member;
    in_file.unsetf(ios::skipws);
    return (in_file);    
}

class Filter_pState : public Filter_State{
public:
    Filter_pState(void) : Filter_State() { }
};

class Filter_cState : public Filter_State{
public:
    Filter_cState(void) : Filter_State() { }
};
#endif
