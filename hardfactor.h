#ifndef _HARD_FACTOR_H_
#define _HARD_FACTOR_H_

#include <algorithm>
#include <list>
#include <map>
#include <string>
#include "integrationcontext.h"
#include "integrationtype.h"

class HardFactorTerm;

class HardFactor {
public:
    virtual const char* get_name() const = 0;
    virtual const size_t get_term_count() const = 0;
    virtual const HardFactorTerm* const* get_terms() const = 0;
};

class HardFactorTerm : public HardFactor {
public:
    HardFactorTerm() : p_this(this) {};
    virtual const IntegrationType* get_type() const = 0;
    virtual void Fs(const IntegrationContext* ictx, double* real, double* imag) const { *real = 0; *imag = 0; }
    virtual void Fn(const IntegrationContext* ictx, double* real, double* imag) const { *real = 0; *imag = 0; }
    virtual void Fd(const IntegrationContext* ictx, double* real, double* imag) const { *real = 0; *imag = 0; }
    const size_t get_term_count() const {
        return 1;
    }
    const HardFactorTerm* const* get_terms() const {
        return &p_this;
    }
private:
    const HardFactorTerm* p_this; // wtf why do I have to do this
};

class HardFactorRegistry {
public:
    void add_hard_factor(const HardFactor* hf, bool manage=false) {
        add_hard_factor(hf->get_name(), hf, manage);
    }
    void add_hard_factor(const char* key, const HardFactor* hf, bool manage=false) {
        using namespace std;
        string keystring(key);
        transform(keystring.begin(), keystring.end(), keystring.begin(), ::tolower);
        hardfactors[keystring] = hf;
        if (manage) {
            hardfactors_to_delete.push_back(hf);
        }
    }
    const HardFactor* get_hard_factor(std::string& key) const {
        using namespace std;
        string keystring(key);
        transform(keystring.begin(), keystring.end(), keystring.begin(), ::tolower);
        map<const string, const HardFactor*>::const_iterator it = hardfactors.find(key);
        if (it == hardfactors.end()) {
            return NULL;
        }
        else {
            return it->second;
        }
    }
protected:
    ~HardFactorRegistry() {
        hardfactors.clear();
        for (std::list<const HardFactor*>::iterator it = hardfactors_to_delete.begin(); it != hardfactors_to_delete.end(); it++) {
            delete (*it);
        }
        hardfactors_to_delete.clear();
    }
private:
    std::map<const std::string, const HardFactor*> hardfactors;
    std::list<const HardFactor*> hardfactors_to_delete;
};

#endif // _HARD_FACTOR_H_
