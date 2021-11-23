#ifndef __DEBUGLVL_H__
#define __DEBUGLVL_H__

class DebugLevel {
   private:
    inline static int _debug_lvl{-1};

   public:
    DebugLevel() = default;
    ~DebugLevel() = default;
    friend void set_debug_lvl(int _lvl);
    friend bool debug_lvl(int _lvl);
};

inline bool debug_lvl(int _lvl) {
    return (DebugLevel::_debug_lvl > _lvl);
}

inline void set_debug_lvl(int _lvl) {
    DebugLevel::_debug_lvl = _lvl;
};

#endif  // __DEBUGLVL_H__