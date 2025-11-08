#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <__msvc_ostream.hpp>
#include <cstdint>

using namespace std;

static const uint8_t    N = 6;             // Size of the board
static const int        INF = 1e9;

struct HistoryEntry {
    uint8_t     player;
    uint8_t     fromCell;
    uint8_t     toCell;
    bool        isJump;
    vector<pair<uint8_t,uint8_t>> flips;
    uint64_t    prevHash;             
};

array<uint8_t, N*N> board;              // 0 - empty, 1 - player 1, 2 - player 2
uint8_t             myPlayer = 1;       
uint8_t             oppPlayer = 2;
uint8_t             pieceCount[3] = {0, 0, 0};
#ifdef DepthM
uint8_t             maxDepth = 6;
#else
uint8_t             maxDepth = 4;
#endif


vector<HistoryEntry> history;

inline bool                     is_inside(int8_t r, int8_t c)   { return r >= 0 && r < N && c >= 0 && c < N; }
inline int8_t                   get_index(int8_t r, int8_t c)   { return r * N + c; }
inline pair<int8_t,int8_t>      get_coordinates(int8_t x)       { return {x / N, x % N}; }

#pragma  region TT
static uint64_t Z_cell[N*N][3];     // 0..2 
static uint64_t Z_side;             // player 2 mask
static uint64_t currentHash = 0;    // curPos hash

static inline uint64_t splitmix64(uint64_t& x){
    uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static void initZobrist(){
    uint64_t seed = 0x12345678ABCDEF01ULL;
    for(uint8_t i=0;i<N * N;++i){
        for(uint8_t v=0; v<3; ++v){
            Z_cell[i][v] = splitmix64(seed);
        }
    }
    Z_side = splitmix64(seed);
}

// Recompute hash if needed
static uint64_t recomputeHash() {
    uint64_t h = 0;
    pieceCount[0] = pieceCount[1] = pieceCount[2] = 0;
    for (uint8_t i = 0; i < N * N; ++i)
    {
        uint8_t v = board[i];
        ++pieceCount[v];
        h ^= Z_cell[i][v];
    }
    currentHash = h;
    return h;
}

inline uint8_t countPieces(uint8_t player) { return pieceCount[player]; }

inline uint8_t countEmpty() { return pieceCount[0]; }


static inline uint64_t posKey(uint8_t player){
    return currentHash ^ (player==2 ? Z_side : 0ULL);
}

enum class TTFlag : uint8_t { EXACT=0, LOWER=1, UPPER=2 };

struct TTEntry {
    uint64_t    key;
    int8_t      depth;   // глубина 
    int32_t     value;   // оценка 
    int8_t      fromCell;
    int8_t      toCell;
    TTFlag      flag;
};

static const size_t TT_SIZE = 1u<<20; 
static std::vector<TTEntry> TT(TT_SIZE, TTEntry{0, -1, 0, -1, -1, TTFlag::EXACT});

static inline TTEntry* ttProbe(uint64_t key){
    return &TT[key & (TT_SIZE-1)];
}
static inline void ttStore(uint64_t key, uint8_t depth, int32_t value, TTFlag flag, uint8_t fromC, uint8_t toC){
    TTEntry* e = ttProbe(key);

    // Store entry with max depth
    if(e->depth <= depth || e->key==0 || e->key==key)
    {
        *e = TTEntry
        {
            .key = key,
            .depth = static_cast<int8_t>(depth),
            .value = value,
            .fromCell = static_cast<int8_t>(fromC),
            .toCell = static_cast<int8_t>(toC),
            .flag = flag
        };
    }
}
#pragma endregion 

// Returns string presentation of cell (Letters - collumns, Digits - rows)
string cellToStr(int8_t idx)
{
    auto [r, c] = get_coordinates(idx);
    string s;
    s += char('a' + c);
    s += char('1' + r);
    return s;
}

// Returns cell index from string presentation
uint8_t strToCell(const string &s)
{
    if (s.size() < 2)
        return -1;
    uint8_t c = s[0] - 'a';
    uint8_t r = s[1] - '1';
    if (!is_inside(r,c))
        return -1;
    return get_index(r,c);
}

// Print board with string presentations of cells
void printBoard(ostream &out)
{
    out << "   a b c d e f\n";
    for (uint8_t r = 0; r < N; ++r)
    {
        out << (r+1) << "  ";
        for (uint8_t c = 0; c < N; ++c)
        {
            out << int(board[get_index(r,c)]);
            if (c+1 < N)
                out << " ";
        }
        out << "\n";
    }
    out << flush;
}

// Initialize board with filled corners
void initBoard()
{
    board.fill(0);

    board[get_index(0,0)] = 1; // a1
    board[get_index(0,N-1)] = 2; // f1
    board[get_index(N-1,0)] = 2; // a6
    board[get_index(N-1,N-1)] = 1; // f6
}

struct Move
{
    int8_t fromCell;
    int8_t toCell;
};

// Collect all possible moves for player
vector<Move> genMoves(uint8_t player)
{
    vector<Move> res;
    res.reserve(128);

    // Find all player's cells
    for (int8_t i = 0; i < N*N; ++i)
    {
        if (board[i] != player)
            continue;
        auto [r, c] = get_coordinates(i);

        // Check all directions
        for (int8_t dr = -2; dr <= 2; ++dr)
        {
            for (int8_t dc = -2; dc <= 2; ++dc)
            {
                if (dr == 0 && dc == 0)
                    continue;
                int8_t nr = r + dr;
                int8_t nc = c + dc;
                if (!is_inside(nr,nc))
                    continue;
                
                int8_t j = get_index(nr,nc);
                if (board[j] != 0)
                    continue;
                
                res.push_back({.fromCell = i, .toCell = j});
            }
        }
    }
    return res;
}

// Apply a move for the given player, updating the board state and recording the move history
void applyMove(uint8_t player, const Move &mv, HistoryEntry &h)
{
    h.player   = player;
    h.fromCell = mv.fromCell;
    h.toCell   = mv.toCell;
    h.flips.clear();
    h.prevHash = currentHash;  

    auto [r1, c1] = get_coordinates(mv.fromCell);
    auto [r2, c2] = get_coordinates(mv.toCell);
    int8_t dist = max(abs(r1-r2), abs(c1-c2));
    h.isJump = (dist == 2);

    if (h.isJump)
    {
        currentHash ^= Z_cell[mv.fromCell][board[mv.fromCell]];
        --pieceCount[board[mv.fromCell]];   
        board[mv.fromCell] = 0;
        ++pieceCount[0];
        currentHash ^= Z_cell[mv.fromCell][board[mv.fromCell]];

        currentHash ^= Z_cell[mv.toCell][board[mv.toCell]];
        --pieceCount[board[mv.toCell]];     
        board[mv.toCell] = player;
        ++pieceCount[player];
        currentHash ^= Z_cell[mv.toCell][board[mv.toCell]];
    }
    else
    {
        currentHash ^= Z_cell[mv.toCell][board[mv.toCell]];
        --pieceCount[board[mv.toCell]];     
        board[mv.toCell] = player;
        ++pieceCount[player];
        currentHash ^= Z_cell[mv.toCell][board[mv.toCell]];
    }

    // flip neighbors
    for (int8_t dr = -1; dr <= 1; ++dr)
    {
        for (int8_t dc = -1; dc <= 1; ++dc)
        {
            if (dr == 0 && dc == 0)
                continue;
            int8_t nr = r2 + dr;
            int8_t nc = c2 + dc;
            if (!is_inside(nr,nc))
                continue;
            int8_t idx = get_index(nr,nc);
            if (board[idx] != 0 && board[idx] != player)
            {
                h.flips.push_back({idx, board[idx]});
                
                currentHash ^= Z_cell[idx][board[idx]];
                --pieceCount[board[idx]];           // opp -> player
                board[idx] = player;
                ++pieceCount[player];
                currentHash ^= Z_cell[idx][board[idx]];
            }
        }
    }
}


// Undo last move
void undoLast()
{
    if (history.empty())
        return;
    HistoryEntry h = history.back();
    history.pop_back();

    // undo flips
    for (auto &p : h.flips) {
        currentHash ^= Z_cell[p.first][board[p.first]];
        --pieceCount[board[p.first]];       
        board[p.first] = p.second;          
        ++pieceCount[board[p.first]];
        currentHash ^= Z_cell[p.first][board[p.first]];
    }

    // undo move
    if (h.isJump)
    {
        currentHash ^= Z_cell[h.toCell][board[h.toCell]];
        --pieceCount[board[h.toCell]];      
        board[h.toCell] = 0;
        ++pieceCount[0];
        currentHash ^= Z_cell[h.toCell][board[h.toCell]];

        currentHash ^= Z_cell[h.fromCell][board[h.fromCell]];
        --pieceCount[board[h.fromCell]];    
        board[h.fromCell] = h.player;
        ++pieceCount[h.player];
        currentHash ^= Z_cell[h.fromCell][board[h.fromCell]];
    }
    else
    {
        currentHash ^= Z_cell[h.toCell][board[h.toCell]];
        --pieceCount[board[h.toCell]];      
        board[h.toCell] = 0;
        ++pieceCount[0];
        currentHash ^= Z_cell[h.toCell][board[h.toCell]];
    }
}

#pragma region Heuristics

static const int8_t D1[8][2] = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,1},{1,-1},{1,0},{1,1}};

int frontierCount(uint8_t player)
{
    int cnt=0;
    for(uint8_t i=0;i<N*N;++i)
    {
        if(board[i]==player)
        {
            uint8_t r=i/N, c=i%N;
            bool front=false;
            for(auto& d : D1)
            {
                int8_t nr=r+d[0], nc=c+d[1];
                if (is_inside(nr,nc) && board[nr*N+nc]==0)
                {
                    front=true;
                    break;
                }
            }
            if(front)
                ++cnt;
        }
    }
    return cnt;
}

// Уязвимые мои фишки: могут быть перевёрнуты за следующий ход соперника
bool oppCanLandNear(int8_t rr, int8_t cc,int opp){
    for(int8_t r=rr-2;r<=rr+2;++r) for(int8_t c=cc-2;c<=cc+2;++c){
        if(!is_inside(r,c))
            continue;
        if(board[r*N+c]!=opp)
            continue;              
        for(int8_t tr=r-2; tr<=r+2; ++tr) for(int8_t tc=c-2; tc<=c+2; ++tc)
        {
            if(!is_inside(tr,tc))
                continue;
            int dist = std::max(abs(tr-r),abs(tc-c));
            if(dist==0 || dist>2)
                continue;
            if(board[tr*N+tc]==0)
            {
                if(std::max(abs(tr-rr),abs(tc-cc))==1)
                    return true;
            }
        }
    }
    return false;
}

int vulnerableCount(int player,int opp){
    int v=0;
    for(int i=0;i<N*N;++i)
    {
        if(board[i]==player)
        {
            int r=i/N, c=i%N;
            if (oppCanLandNear(r,c,opp))
                ++v;
        }
    }
    return v;
}

int liberties(int player){
    int sum=0;
    for(int i=0;i<N*N;++i) if(board[i]==player){
        int r=i/N, c=i%N;
        for(auto& d : D1){
            int nr=r+d[0], nc=c+d[1];
            if(is_inside(nr,nc) && board[nr*N+nc]==0) ++sum;
        }
    }
    return sum;
}

#pragma endregion 

int evaluateSimple()
{
    uint8_t myCnt = countPieces(myPlayer);
    uint8_t oppCnt = countPieces(oppPlayer);
   
#ifdef Weightx10
    int16_t count = (myCnt - oppCnt) * 10;
#elif defined Weightx5
    int16_t count = (myCnt - oppCnt) * 5;
#else
    int16_t count = (myCnt - oppCnt);
#endif

    int frontier = 0;
    int vul = 0;
    int liber = 0;
    

#ifdef Frontier
    int f_my   = frontierCount(myPlayer);
    int f_opp  = frontierCount(oppPlayer);
    frontier = f_opp - f_my;
#endif

#ifdef Vuln
    int vul_my = vulnerableCount(myPlayer, oppPlayer);
    int vul_op = vulnerableCount(oppPlayer, myPlayer);
    vul = vul_op - vul_my;
#endif

#ifdef Liber
    int lib_my = liberties(myPlayer);
    int lib_op = liberties(oppPlayer);
    liber = lib_my - lib_op;
#endif

#ifdef DepthM
    auto movesCount = 0;
#else
    int myMoves = (genMoves(myPlayer).size());
    int opMoves = (genMoves(oppPlayer).size());
    int movesCount = myMoves - opMoves;
#endif

    uint8_t empties = countEmpty();
    double phase = 1.0 - static_cast<double>(empties) / (N*N);

    //int16_t score = count;
    double score =
        1.0 * count
      + 0.8 * movesCount
      + 0.6 * liber
      + 0.8 * frontier          
      + 0.9 * vul       
      ;
    
#ifdef Phase
    score = score*(1.0 - 0.7*phase) + (count)*0.7*phase;
#endif
    
    return std::round(score);
}



static int negamax(uint8_t depth, int alpha, int beta, uint8_t player)
{
    const int alphaOrig = alpha;
    uint64_t key = posKey(player);

    // TT probe
    TTEntry* te = ttProbe(key);
    if (te->key == key && te->depth >= depth) {
        int v = te->value;
        if (te->flag == TTFlag::EXACT)
            return v;

        if (te->flag == TTFlag::LOWER)
            alpha = std::max(alpha, v);
        else if (te->flag == TTFlag::UPPER)
            beta  = std::min(beta,  v);

        if (alpha >= beta)
            return v;
    }

    auto moves = genMoves(player);

    // End Of Game
    if (moves.empty())
    {
        int myCnt  = countPieces(myPlayer);
        int oppCnt = countPieces(oppPlayer);
        int winner = (myCnt > oppCnt) ? +1 : (myCnt < oppCnt ? -1 : 0);
        int val = winner * 100000;
        if (player != myPlayer) 
            val = -val;
        ttStore(key, depth, val, TTFlag::EXACT, -1, -1);
        return val;
    }

    // Leaf
    if (depth == 0)
    {
        int val = evaluateSimple();
        if (player != myPlayer)
            val = -val;
        ttStore(key, depth, val, TTFlag::EXACT, -1, -1);
        return val;
    }

    // Check if best move was found
    if (te->key == key && te->fromCell >= 0 && te->toCell >= 0)
    {
        auto it = std::find_if(moves.begin(), moves.end(), [&](const Move& m)
        {
            return m.fromCell==te->fromCell && m.toCell==te->toCell;
        });
        if (it!=moves.end())
            std::rotate(moves.begin(), it, it+1);
    }

    int best = -INF;
    Move bestMove{-1,-1};

    for (auto &mv : moves)
    {
        // Apply move
        HistoryEntry h;
        applyMove(player, mv, h);
        history.push_back(h);

        // Go down
        int val = -negamax(depth-1, -beta, -alpha, (player == 1) ? 2 : 1);

        // Undo move
        undoLast();

        // Is the move best?
        if (val > best)
        {
            best = val;
            bestMove = mv;
        }
        
        alpha = std::max(alpha, val);
        if (alpha >= beta)
            break;
    }

    // TT store
    TTFlag flag = TTFlag::EXACT;
    if (best <= alphaOrig)
        flag = TTFlag::UPPER;
    else if (best >= beta)
        flag = TTFlag::LOWER;
    ttStore(key, depth, best, flag, bestMove.fromCell, bestMove.toCell);

    return best;
}


static bool chooseBestMove(Move &bestMove)
{
    auto rootMoves = genMoves(myPlayer);
    if (rootMoves.empty())
        return false;

    uint64_t k = posKey(myPlayer);
    TTEntry* e = ttProbe(k);
    if (e->key==k && e->fromCell>=0)
    {
        auto it = std::find_if(rootMoves.begin(), rootMoves.end(),
            [&](const Move& x){ return x.fromCell==e->fromCell && x.toCell==e->toCell; });
        if (it != rootMoves.end())
            std::rotate(rootMoves.begin(), it, it+1);
    }

    // Best - first
    auto promoteBestFirst = [&](const Move& m){
        auto it = std::find_if(rootMoves.begin(), rootMoves.end(),
                               [&](const Move& x){ return x.fromCell==m.fromCell && x.toCell==m.toCell; });
        if (it!=rootMoves.end()) {
            std::rotate(rootMoves.begin(), it, it+1); // found move in begin
        }
    };

    Move globalBest   = rootMoves[0];
    
    for (int depth = 1; depth <= maxDepth; ++depth)
    {
        // first - best
        if (depth > 1)
            promoteBestFirst(globalBest);

        int alpha = -INF;
        int beta  =  INF;

        int bestValThisDepth = -INF;
        Move bestMoveThisDepth = rootMoves[0];

        auto searchWithWindow = [&](int a, int b, int& outVal, Move& outMove) -> bool
        {
            int localBest = -INF;
            Move localBestMove = rootMoves[0];

            for (auto &mv : rootMoves)
            {
                HistoryEntry h;
                applyMove(myPlayer, mv, h);
                history.push_back(h);

                int val = -negamax(depth-1, -b, -a, oppPlayer);

                undoLast();

                if (val > localBest) {
                    localBest = val;
                    localBestMove = mv;
                }
                a = std::max(a, val);
                if (a >= b)
                    break; 
            }
            outVal = localBest;
            outMove = localBestMove;
            return true;
        };

        int score = bestValThisDepth;
        searchWithWindow(alpha, beta, score, bestMoveThisDepth);

        globalBest    = bestMoveThisDepth;
    }

    bestMove = globalBest;
    return true;
}


// Check if move is legal for given player
bool isLegalMove(int player, const Move& mv) {
    if (mv.fromCell < 0 || mv.fromCell >= N*N)
        return false;
    if (mv.toCell   < 0 || mv.toCell   >= N*N)
        return false;

    if (board[mv.fromCell] != player)
        return false; 
    if (board[mv.toCell]   != 0)
        return false; 

    auto [r1, c1] = get_coordinates(mv.fromCell);
    auto [r2, c2] = get_coordinates(mv.toCell);
    if (!is_inside(r1,c1) || !is_inside(r2,c2))
        return false;

    int dist = max(abs(r1 - r2), abs(c1 - c2));
    if (dist != 1 && dist != 2)
        return false;     

    return true;
}

void finishAndExit()
{
    int myCnt = countPieces(myPlayer);
    int oppCnt = countPieces(oppPlayer);

    if (myCnt > oppCnt)
        exit(0);
    else if (myCnt < oppCnt)
        exit(3);
    else
        exit(4);
}

void printMyMoveAndBoard(const Move &mv)
{
    string sFrom = cellToStr(mv.fromCell);
    string sTo   = cellToStr(mv.toCell);

    // botctl
    cerr << sFrom << " " << sTo << "\n";

    // console
    cout << "Bot move: " << sFrom << " " << sTo << "\n";
    printBoard(cout);
}

bool isEmptyCells()
{
    for (int i = 0; i < N*N; ++i)
        if (board[i] == 0)
            return true;
    return false;
}

bool oppHasMove(int player) {
    // player — это соперник (например oppPlayer)
    for (int i = 0; i < N*N; ++i) {
        if (board[i] != player) continue;
        auto [r, c] = get_coordinates(i);
        for (int dr = -2; dr <= 2; ++dr) {
            for (int dc = -2; dc <= 2; ++dc) {
                if (dr == 0 && dc == 0) continue;
                int nr = r + dr;
                int nc = c + dc;
                if (!is_inside(nr, nc)) continue;
                int j = get_index(nr, nc);
                if (board[j] == 0) {
                    int dist = max(abs(dr), abs(dc));
                    if (dist == 1 || dist == 2) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}


int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    initBoard();
    initZobrist();
    recomputeHash();

    for (int i=0;i<argc;++i) std::cout << "argv["<<i<<"]="<<argv[i]<<"\n";

    // player arg - last
    if (argc >= 2)
    {
        string last = argv[argc-1];
        if (last == "0" || last == "1")
        {
            int clr = last[0] - '0'; 
            myPlayer = clr + 1;      
        }
        // depth
        for (int i = 1; i+1 < argc; ++i)
        {
            string a = argv[i];
            if (a == "-d" && i+2 <= argc)
            {
                maxDepth = max(1, atoi(argv[i+1]));
            }
        }
    }
    else
    {
        myPlayer = 1;
    }
    oppPlayer = (myPlayer == 1) ? 2 : 1;

    // bot - first
    if (myPlayer == 1)
    {
        Move mv;
        if (!chooseBestMove(mv))
        {
            finishAndExit();
        }
        HistoryEntry h;
        applyMove(myPlayer, mv, h);
        history.push_back(h);
        printMyMoveAndBoard(mv);
    }

    // game cycle
    while (true)
    {
        string tok;
        if (!(cin >> tok))
        {
            finishAndExit();
        }

        if (tok == "u1") {
            // undo 1 player move and 1 bot move
            if (!history.empty())
                undoLast();
            if (!history.empty())
                undoLast();
            
            cout << "Undo performed.\n";
            printBoard(cout);
            continue;
        }

        // Apply opponent's move
        string tok2;
        if (!(cin >> tok2))
        {
            finishAndExit();
        }
        int8_t from = strToCell(tok);
        int8_t to   = strToCell(tok2);
        if (from == -1 || to == -1)
        {
            cout << "Illegal move: " << tok << " " << tok2 << "\n";
            printBoard(cout);
            continue;
        }
        
        Move mv{from, to};
        if (!isLegalMove(oppPlayer, mv)) {
            cout << "Illegal move: " << tok << " " << tok2 << "\n";
            printBoard(cout);
            continue; 
        }
        HistoryEntry h;
        applyMove(oppPlayer, mv, h);
        history.push_back(h);
            
        cout << "Opponent (player) move: " << tok << " " << tok2 << "\n";
        printBoard(cout);

        // Bot's move
        Move mymv;
        if (!chooseBestMove(mymv))
        {
            finishAndExit();
        }

        {
            HistoryEntry h;
            applyMove(myPlayer, mymv, h);
            history.push_back(h);
            printMyMoveAndBoard(mymv);
        }

        if (!isEmptyCells())
            finishAndExit();
        if (!oppHasMove(oppPlayer))
            finishAndExit();
    }

    return 0;
}
