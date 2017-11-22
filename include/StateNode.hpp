#ifndef STATE_NODE_HPP
#define STATE_NODE_HPP
#include <solution.h>

template<typename T>
class StateNode {
    private:
        T obj;
        Job *job;
        StateNode *pred;

        inline void getData (T &o, Job* &j, StateNode* &p) const {
            o = obj;
            j = job;
            p = pred;
        }

    public:
        /*--------------------------
         | standard query functions
         *-------------------------*/

        inline Job* getJob() const {return job;}
        inline T getObj() const {return obj;}
        inline StateNode *getPred() const {return pred;}


        /*--------------------------------
         | copy another StateNode to this one
         *-------------------------------*/

        inline void copy (StateNode *s) {s->getData(obj, job, pred);}

        inline void copy (StateNode &s) {s.getData(obj, job, pred);}

        inline void init (T _obj, Job  *j, StateNode *p) {
            job = j;
            obj = _obj;
            pred   = p;
        }

        inline void extend (Job *_job,  double cost, StateNode *source) {
            pred   = source;                        //source will be predecessor
            obj = source->getObj() + cost; //length takes into account the arc
            job = _job;
        }

        void output (FILE *file, const char *before=NULL, const char *after=NULL) {
            if (before) fprintf (file, before);
            fprintf (file, "(%d", job->job);
            fprintf (file, ",%.1f)", obj);
            if (after) fprintf (file, after);
        }

        void outputPath(FILE *file, int k) {
            if (getPred() && k>1) {
                getPred()->outputPath(file, k-1);
                fprintf (file, ":");
            }
            fprintf (file, "%d", job->job);
        }

        inline bool sameEnding (StateNode *s, StateNode *t, int k) const {
            while (s && t) {
                if (s->getJob()!=t->getJob()) return false;
                if (--k == 0) return true;
                s = s->getPred();
                t = t->getPred();
            }
            return (s==t);
        }

        inline void markForbidden (bool *forbidden, bool value, int k) {
            StateNode *s = this;
            do { //assuming c is at least 1!
                forbidden[s->getJob()] = value;
                s = s->getPred();
                k--;
            } while (s && k);
        }

        inline StateNode *getSemiDominator (StateNode *s, int k) {
            if (sameEnding(s,this,k)) return (s->getObj()>obj ? s : this);
            return NULL;
        }


        inline int getPath (Job *path, int k) {
            int len = 0;
            StateNode *s = this;
            do {
                path[len++] = s->getJob();
                k--;
                s = s->getPred();
            } while (k && s);
            return len;
        }
};



#endif // STATE_NODE_HPP
