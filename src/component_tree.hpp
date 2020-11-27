#pragma once
#include <opencv2/core.hpp>
#include <iostream>

// namespace morpho
//{
// namespace internal
//{
// template<bool parallel = false>
// struct maxtree_flood_algorithm
//{
//    static constexpr std::size_t nlevels = 1ul << value_traits<uint>::quant;
//    static constexpr bool use_dejavu     = true;
//
//    static constexpr uint UNINITIALIZED = std::numeric_limits<uint>::max();
//    static constexpr uint INQUEUE       = 0;
//
//    uint flood(uint level)
//    {
//        mln_precondition(m_has_repr[level]);
//        mln_precondition(not m_q.empty(level));
//        uint r = m_repr[level];
//        mln_precondition(m_h(m_ima[r]) == level);
//
//        // std::cout << "Start of flood: " << (int)level << "(repr:" << r << ")" << std::endl;
//        // io::imprint(m_parent);
//
//        while (!m_q.empty(level))
//        {
//            uint p = m_q.pop_at_level(level);
//            if (!parallel and p != r)
//            {
//                *(--m_out) = p;
//            }
//            m_parent[p] = r;
//
//            mln_foreach(auto k, m_nbh_delta_indexes)
//            {
//                auto q = p + k;
//
//                bool processed;
//                if (use_dejavu)
//                    processed = m_deja_vu[q];
//                else if (!parallel) // no bound checking
//                    processed = (m_parent[q] != UNINITIALIZED);
//                else
//                    processed = q < m_first_index or m_last_index <= q or (m_parent[q] != UNINITIALIZED);
//
//                if (!processed)
//                {
//                    uint newlevel = m_h(m_ima[q]);
//                    if (!m_has_repr[newlevel])
//                    {
//                        m_repr[newlevel]     = q;
//                        m_has_repr[newlevel] = true;
//                    }
//                    // std::cout << "++ push: " << q << " @ " << (int)newlevel << std::endl;
//                    m_q.push_at_level(q, newlevel);
//                    if (use_dejavu)
//                        m_deja_vu[q] = true;
//                    else
//                        m_parent[q] = INQUEUE;
//
//                    if (level < newlevel)
//                        do
//                        {
//                            newlevel = flood(newlevel);
//                        } while (level < newlevel);
//                }
//            }
//        }
//
//        // Attach to parent
//        if (!parallel)
//        {
//            *(--m_out) = r;
//        }
//        m_has_repr[level] = false;
//        while (level > value_traits<uint>::min())
//        {
//            --level;
//            if (m_has_repr[level])
//            {
//                m_parent[r] = m_repr[level];
//                break;
//            }
//        }
//
//        return level;
//    }
//
//    maxtree_flood_algorithm(const cv::Mat& ima, cv::Mat& parent, uint* Send = NULL)
//        : m_ima(ima),
//          m_parent(parent),
//          m_has_repr {
//              false,
//          },
//          m_out(Send)
//    {
//        if (use_dejavu)
//        {
//            resize(m_deja_vu, ima).init(false);
//            extension::fill(m_deja_vu, true);
//        }
//
//        m_first_index = m_ima.index_of_point(m_ima.domain().pmin);
//        m_last_index  = m_ima.index_of_point(m_ima.domain().pmax) - ima.index_strides()[0];
//
//        // Get min element and reserve queue
//        size_t pmin = ima.index_of_point(ima.domain().pmin);
//        uint vmin   = value_traits<uint>::max();
//        {
//            std::vector<std::size_t> h(nlevels, 0);
//            // uint h[nlevels] = {0,};
//
//            mln_pixter(px, ima);
//            mln_forall(px)
//            {
//                uint l = m_h(px->val());
//                ++h[l];
//                if (l < vmin)
//                {
//                    vmin = l;
//                    pmin = px->index();
//                }
//            }
//
//            m_q.init(&h[0]);
//        }
//
//        // Start flooding
//        // std::cout << (int) vmin << "@" << pmin << std::endl;
//        m_q.push_at_level(pmin, vmin);
//        m_repr[vmin]     = pmin;
//        m_has_repr[vmin] = true;
//        if (use_dejavu)
//            m_deja_vu[pmin] = true;
//        else
//            m_parent[pmin] = INQUEUE;
//        flood(vmin);
//    }
//
//    static void run(const cv::Mat& ima, cv::Mat& parent, uint* Send)
//    {
//        maxtree_flood_algorithm x(ima, parent, nbh, cmp, Send);
//        if (!parallel)
//            assert((x.m_out + ima.domain().size()) == Send);
//        (void)x;
//    }
//
// private:
//    typedef mln_ch_value(I, bool) J;
//
//    StrictWeakOrdering m_cmp;
//    indexer<V, StrictWeakOrdering> m_h;
//
//    const I& m_ima;
//    mln_ch_value(I, uint) & m_parent;
//    mln_ch_value(I, bool) m_deja_vu;
//    std::array<typename I::difference_type, Neighborhood::static_size> m_nbh_delta_indexes;
//
//    bounded_hqueue<uint, nlevels> m_q;
//    bool m_has_repr[nlevels];
//    uint m_repr[nlevels];
//    uint* m_out;
//    uint m_first_index;
//    uint m_last_index;
//};
//
//} // namespace internal
// template<bool parallel = false>
// struct MaxTreeAlgorithmHQ
//{
//    static constexpr const uint UNINITIALIZED = std::numeric_limits<uint>::max();
//    static constexpr const uint INQUEUE       = 0;
//    static constexpr const bool use_dejavu    = false;
//
//    MaxTreeAlgorithmHQ(const cv::Mat& ima): m_ima(ima), m_has_previous(false)
//    {
//        m_parent = cv::Mat(ima.size(), CV_32S);
//
//        if (!use_dejavu)
//        {
//            m_parent.setTo(cv::Scalar(UNINITIALIZED));
//            // resize(m_parent, ima).init((uint)UNINITIALIZED);
//            // extension::fill(m_parent, (uint)INQUEUE);
//        }
//
//        m_nsplit = 0;
//        if (!parallel)
//        {
//            m_S.resize(ima.cols * ima.rows);
//        }
//    }
//
//
//    MaxTreeAlgorithmHQ(MaxTreeAlgorithmHQ& other): m_ima(other.m_ima), m_parent(other.m_parent), m_has_previous(false)
//    {
//        m_nsplit = 0;
//    }
//
//
//    void run(cv::Rect domain)
//    {
//        cv::Mat ima    = m_ima(domain);
//        cv::Mat parent = m_parent(domain);
//
//        if (parallel)
//            internal::maxtree_flood_algorithm<parallel>::run(ima, parent, NULL);
//        else
//            internal::maxtree_flood_algorithm<parallel>::run(ima, parent, &m_S[0] + domain.area());
//
//        if (m_has_previous)
//        {
//            // this->join(*this, false);
//            // m_current_domain.join(domain);
//            // m_nsplit += 1;
//        }
//        else
//        {
//            m_current_domain = domain;
//            m_has_previous   = true;
//        }
//    }
//
//    // void join(MaxTreeAlgorithmHQ& other, bool joindomain = true)
//    //{
//    //    mln_precondition(m_has_previous);
//
//    //    merge_tree(m_ima, m_parent, this->m_current_domain, m_cmp);
//    //    if (joindomain)
//    //    {
//    //        m_current_domain.join(other.m_current_domain);
//    //        m_nsplit += other.m_nsplit + 1;
//    //    }
//    //}
//
//
// public:
//    const cv::Mat& m_ima;
//    cv::Rect m_current_domain;
//    cv::Mat m_parent;
//    bool m_has_previous;
//    unsigned m_nsplit;
//    std::vector<uint> m_S;
//};
//
//
//// namespace parallel
////{
//// template<typename V, typename Neighborhood, typename StrictWeakOrdering = std::less<V>>
//// std::pair<image2d<typename image2d<V>::uint>, std::vector<typename image2d<V>::uint>>
////    maxtree_hqueue(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp = StrictWeakOrdering())
////{
////    typedef typename image2d<V>::uint uint;
////    MaxTreeAlgorithmHQ<V, Neighborhood, StrictWeakOrdering, true> algo(ima, nbh, cmp);
////    int nmaxsplit = tbb::task_scheduler_init::default_num_threads() * 4;
////    int grain     = std::max(ima.nrows() / nmaxsplit, 1u);
////    std::cout << "Grain: " << grain << std::endl;
////    tbb::parallel_reduce(grain_box2d(ima.domain(), grain), algo, tbb::auto_partitioner());
////
////    std::cout << "Number of split: " << algo.m_nsplit << std::endl;
////    image2d<uint>& parent = algo.m_parent;
////    std::vector<uint> S(ima.domain().size());
////    canonize(ima, parent, &S[0]);
////
////    // MaxtreeCanonizationAlgorithm<V> canonizer(ima, parent);
////    // tbb::parallel_for(grain_box2d(ima.domain(), grain), canonizer, tbb::auto_partitioner());
////
////    return std::make_pair(std::move(parent), std::move(S));
////}
////} // namespace parallel
//
//
// namespace serial
//{
// void maxtree_hqueue(const cv::Mat& ima, cv::Mat& parent, std::vector<uint>& S)
//{
//    MaxTreeAlgorithmHQ<false> algo(ima);
//    algo.run();
//    std::cout << "Number of split: " << algo.m_nsplit << std::endl;
//
//    parent = algo.m_parent;
//    S      = algo.m_S;
//}
//
//} // namespace serial
//
//
//} // namespace morpho

using std::string;
using std::vector;


using namespace cv;
struct tree
{
    tree() {}
    struct Params
    {
        Params(int _delta             = 5,
               int _min_area          = 60,
               int _max_area          = 14400,
               double _max_variation  = 0.25,
               double _min_diversity  = .2,
               int _max_evolution     = 200,
               double _area_threshold = 1.01,
               double _min_margin     = 0.003,
               int _edge_blur_size    = 5)
        {
            delta         = _delta;
            minArea       = _min_area;
            maxArea       = _max_area;
            maxVariation  = _max_variation;
            minDiversity  = _min_diversity;
            maxEvolution  = _max_evolution;
            areaThreshold = _area_threshold;
            minMargin     = _min_margin;
            edgeBlurSize  = _edge_blur_size;
            pass2Only     = false;
        }

        int delta;
        int minArea;
        int maxArea;
        double maxVariation;
        double minDiversity;
        bool pass2Only;

        int maxEvolution;
        double areaThreshold;
        double minMargin;
        int edgeBlurSize;
    };
    enum
    {
        DIR_SHIFT = 29,
        NEXT_MASK = ((1 << DIR_SHIFT) - 1)
    };

    struct Pixel
    {
        Pixel(): val(0) {}
        Pixel(int _val): val(_val) {}

        int getGray(const Pixel* ptr0, const uchar* imgptr0, int mask) const
        {
            return imgptr0[this - ptr0] ^ mask;
        }
        int getNext() const
        {
            return (val & NEXT_MASK);
        }
        void setNext(int next)
        {
            val = (val & ~NEXT_MASK) | next;
        }

        int getDir() const
        {
            return (int)((unsigned)val >> DIR_SHIFT);
        }
        void setDir(int dir)
        {
            val = (val & NEXT_MASK) | (dir << DIR_SHIFT);
        }
        bool isVisited() const
        {
            return (val & ~NEXT_MASK) != 0;
        }

        int val;
    };
    typedef int PPixel;

    struct WParams
    {
        WParams() {};
        Params p;
        vector<vector<cv::Point>>* msers;
        vector<cv::Rect>* bboxvec;
        Pixel* pix0;
        int step;
    };
    // the history of region grown
    struct CompHistory
    {
        CompHistory()
        {
            parent_ = child_ = next_ = 0;
            val = size = 0;
            var        = -1.f;
            head       = 0;
            checked    = false;
        }
        void updateTree(WParams& wp, CompHistory** _h0, CompHistory** _h1, bool final)
        {
            if (var >= 0.f)
                return;
            int delta = wp.p.delta;

            CompHistory *h0_ = 0, *h1_ = 0;
            CompHistory* c = child_;
            // printf("sz: %d\n", size);
            if (size >= wp.p.minArea)
            {
                for (; c != 0; c = c->next_)
                {
                    if (c->var < 0.f)
                        c->updateTree(wp, c == child_ ? &h0_ : 0, c == child_ ? &h1_ : 0, final);
                    if (c->var < 0.f)
                        return;
                }
            }

            // find h0 and h1 such that:
            //    h0->val >= h->val - delta and (h0->parent == 0 or h0->parent->val < h->val - delta)
            //    h1->val <= h->val + delta and (h1->child == 0 or h1->child->val < h->val + delta)
            // then we will adjust h0 and h1 as h moves towards latest
            CompHistory *h0 = this, *h1 = h1_ && h1_->size > size ? h1_ : this;
            if (h0_)
            {
                for (h0 = h0_; h0 != this && h0->val < val - delta; h0 = h0->parent_)
                    ;
            }
            else
            {
                for (; h0->child_ && h0->child_->val >= val - delta; h0 = h0->child_)
                    ;
            }

            for (; h1->parent_ && h1->parent_->val <= val + delta; h1 = h1->parent_)
                ;

            if (_h0)
                *_h0 = h0;
            if (_h1)
                *_h1 = h1;

            // when we do not well-defined ER(h->val + delta), we stop
            // the process of computing variances unless we are at the final step
            if (!final && !h1->parent_ && h1->val < val + delta)
                return;

            // printf("h1 size: %d, h0 size: %d \n", h1->size, h0->size);
            var = (float)(h1->size - h0->size) / size;
            c   = child_;
            for (; c != 0; c = c->next_)
                c->checkAndCapture(wp);
            if (final && !parent_)
                checkAndCapture(wp);
        }

        void checkAndCapture(WParams& wp)
        {
            if (checked)
                return;
            checked = true;
            if (size < wp.p.minArea || size > wp.p.maxArea || var < 0.f || var > wp.p.maxVariation)
                return;
            if (child_)
            {
                CompHistory* c = child_;
                for (; c != 0; c = c->next_)
                {
                    if (c->var >= 0.f && var > c->var)
                        return;
                }
            }
            if (var > 0.f && parent_ && parent_->var >= 0.f && var >= parent_->var)
                return;
            int xmin = INT_MAX, ymin = INT_MAX, xmax = INT_MIN, ymax = INT_MIN, j = 0;
            wp.msers->push_back(vector<cv::Point>());
            vector<cv::Point>& region = wp.msers->back();
            region.resize(size);
            const Pixel* pix0 = wp.pix0;
            int step          = wp.step;

            for (PPixel pix = head; j < size; j++, pix = pix0[pix].getNext())
            {
                int y = pix / step;
                int x = pix - y * step;

                xmin = std::min(xmin, x);
                xmax = std::max(xmax, x);
                ymin = std::min(ymin, y);
                ymax = std::max(ymax, y);

                region[j] = Point(x, y);
            }

            wp.bboxvec->push_back(Rect(xmin, ymin, xmax - xmin + 1, ymax - ymin + 1));
        }

        CompHistory* child_;
        CompHistory* parent_;
        CompHistory* next_;
        int val;
        int size;
        float var;
        PPixel head;
        bool checked;
    };

    struct ConnectedComp
    {
        ConnectedComp()
        {
            init(0);
        }

        void init(int gray)
        {
            head = tail = 0;
            history     = 0;
            size        = 0;
            gray_level  = gray;
        }

        // add history chunk to a connected component
        void growHistory(CompHistory*& hptr, WParams& wp, int new_gray_level, bool final)
        {
            if (new_gray_level < gray_level)
                new_gray_level = gray_level;

            CompHistory* h;
            if (history && history->val == gray_level)
            {
                h = history;
            }
            else
            {
                h          = hptr++;
                h->parent_ = 0;
                h->child_  = history;
                h->next_   = 0;

                if (history)
                {
                    history->parent_ = h;
                }
            }
            CV_Assert(h != NULL);
            h->val     = gray_level;
            h->size    = size;
            h->head    = head;
            h->var     = FLT_MAX;
            h->checked = true;
            if (h->size >= wp.p.minArea)
            {
                h->var     = -1.f;
                h->checked = false;
            }

            gray_level = new_gray_level;
            history    = h;
            if (history && history->val != gray_level)
            {
                history->updateTree(wp, 0, 0, final);
            }
        }

        // merging two connected components
        void merge(ConnectedComp* comp1, ConnectedComp* comp2, CompHistory*& hptr, WParams& wp)
        {
            if (comp1->gray_level < comp2->gray_level)
                std::swap(comp1, comp2);

            gray_level = comp1->gray_level;
            comp1->growHistory(hptr, wp, gray_level, false);
            comp2->growHistory(hptr, wp, gray_level, false);

            if (comp1->size == 0)
            {
                head = comp2->head;
                tail = comp2->tail;
            }
            else
            {
                head = comp1->head;
                wp.pix0[comp1->tail].setNext(comp2->head);
                tail = comp2->tail;
            }

            size    = comp1->size + comp2->size;
            history = comp1->history;

            CompHistory* h1 = history->child_;
            CompHistory* h2 = comp2->history;
            // the child_'s size should be the large one
            if (h1 && h1->size > h2->size)
            {
                // add h2 as a child only if its size is large enough
                if (h2->size >= wp.p.minArea)
                {
                    h2->next_   = h1->next_;
                    h1->next_   = h2;
                    h2->parent_ = history;
                }
            }
            else
            {
                history->child_ = h2;
                h2->parent_     = history;
                // reserve h1 as a child only if its size is large enough
                if (h1 && h1->size >= wp.p.minArea)
                {
                    h2->next_ = h1;
                }
            }
        }

        PPixel head;
        PPixel tail;
        CompHistory* history;
        int gray_level;
        int size;
    };
    void preprocess1(const Mat& img, int* level_size)
    {
        memset(level_size, 0, 256 * sizeof(level_size[0]));

        int i, j, cols = img.cols, rows = img.rows;
        int step = cols;
        pixbuf.resize(step * rows);
        heapbuf.resize(cols * rows + 256);
        histbuf.resize(cols * rows);
        Pixel borderpix;
        borderpix.setDir(5);

        for (j = 0; j < step; j++)
        {
            pixbuf[j] = pixbuf[j + (rows - 1) * step] = borderpix;
        }

        for (i = 1; i < rows - 1; i++)
        {
            const uchar* imgptr = img.ptr(i);
            Pixel* pptr         = &pixbuf[i * step];
            pptr[0] = pptr[cols - 1] = borderpix;
            for (j = 1; j < cols - 1; j++)
            {
                int val = imgptr[j];
                level_size[val]++;
                pptr[j].val = 0;
            }
        }
    }

    void preprocess2(const Mat& img, int* level_size)
    {
        int i;

        for (i = 0; i < 128; i++)
            std::swap(level_size[i], level_size[255 - i]);

        if (!params.pass2Only)
        {
            int j, cols = img.cols, rows = img.rows;
            int step = cols;
            for (i = 1; i < rows - 1; i++)
            {
                Pixel* pptr = &pixbuf[i * step];
                for (j = 1; j < cols - 1; j++)
                {
                    pptr[j].val = 0;
                }
            }
        }
    }

    void pass(const Mat& img,
              vector<vector<Point>>& msers,
              vector<Rect>& bboxvec,
              Size size,
              const int* level_size,
              int mask)
    {
        CompHistory* histptr = &histbuf[0];
        int step             = size.width;
        Pixel *ptr0 = &pixbuf[0], *ptr = &ptr0[step + 1];
        const uchar* imgptr0 = img.ptr();
        Pixel** heap[256];
        ConnectedComp comp[257];
        ConnectedComp* comptr = &comp[0];
        WParams wp;
        wp.p       = params;
        wp.msers   = &msers;
        wp.bboxvec = &bboxvec;
        wp.pix0    = ptr0;
        wp.step    = step;

        heap[0]    = &heapbuf[0];
        heap[0][0] = 0;

        for (int i = 1; i < 256; i++)
        {
            heap[i]    = heap[i - 1] + level_size[i - 1] + 1;
            heap[i][0] = 0;
        }

        comptr->gray_level = 256;
        comptr++;
        comptr->gray_level = ptr->getGray(ptr0, imgptr0, mask);
        ptr->setDir(1);
        int dir[] = {0, 1, step, -1, -step};
        for (;;)
        {
            int curr_gray = ptr->getGray(ptr0, imgptr0, mask);
            int nbr_idx   = ptr->getDir();
            // take tour of all the 4 directions
            for (; nbr_idx <= 4; nbr_idx++)
            {
                // get the neighbor
                Pixel* ptr_nbr = ptr + dir[nbr_idx];
                if (!ptr_nbr->isVisited())
                {
                    // set dir=1, next=0
                    ptr_nbr->val = 1 << DIR_SHIFT;
                    int nbr_gray = ptr_nbr->getGray(ptr0, imgptr0, mask);
                    if (nbr_gray < curr_gray)
                    {
                        // when the value of neighbor smaller than current
                        // push current to boundary heap and make the neighbor to be the current one
                        // create an empty comp
                        *(++heap[curr_gray]) = ptr;
                        ptr->val             = (nbr_idx + 1) << DIR_SHIFT;
                        ptr                  = ptr_nbr;
                        comptr++;
                        comptr->init(nbr_gray);
                        curr_gray = nbr_gray;
                        nbr_idx   = 0;
                        continue;
                    }
                    // otherwise, push the neighbor to boundary heap
                    *(++heap[nbr_gray]) = ptr_nbr;
                }
            }

            // set dir = nbr_idx, next = 0
            ptr->val   = nbr_idx << DIR_SHIFT;
            int ptrofs = (int)(ptr - ptr0);
            CV_Assert(ptrofs != 0);

            // add a pixel to the pixel list
            if (comptr->tail)
                ptr0[comptr->tail].setNext(ptrofs);
            else
                comptr->head = ptrofs;
            comptr->tail = ptrofs;
            comptr->size++;
            // get the next pixel from boundary heap
            if (*heap[curr_gray])
            {
                ptr = *heap[curr_gray];
                heap[curr_gray]--;
            }
            else
            {
                for (curr_gray++; curr_gray < 256; curr_gray++)
                {
                    if (*heap[curr_gray])
                        break;
                }
                if (curr_gray >= 256)
                    break;

                ptr = *heap[curr_gray];
                heap[curr_gray]--;

                if (curr_gray < comptr[-1].gray_level)
                {
                    comptr->growHistory(histptr, wp, curr_gray, false);
                    CV_DbgAssert(comptr->size == comptr->history->size);
                }
                else
                {
                    // there must one pixel with the second component's gray level in the heap,
                    // so curr_gray is not large than the second component's gray level
                    comptr--;
                    CV_DbgAssert(curr_gray == comptr->gray_level);
                    comptr->merge(comptr, comptr + 1, histptr, wp);
                    CV_DbgAssert(curr_gray == comptr->gray_level);
                }
            }
        }

        for (; comptr->gray_level != 256; comptr--)
        {
            comptr->growHistory(histptr, wp, 256, true);
        }
    }

    vector<Pixel> pixbuf;
    vector<Pixel*> heapbuf;
    vector<CompHistory> histbuf;

    Params params;
};