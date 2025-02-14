
#include "labellist.h"

static
int deleteChildren2(
        std::vector<LabelList*>&    activeLists,
        std::vector<LabelList*>&    oldLists,
        LabelNode*                  node,
        LabelNode**                 next_node
){
    LabelNode* child_node;
    if(node->is_propagated_)
    {
        int deleted_labels = 1;
        child_node = node->child_;
        while (child_node != nullptr)
        {
            deleted_labels += deleteChildren2(activeLists, oldLists, child_node, next_node);
            child_node = node->child_;
        }
        *next_node = node->next_;
        oldLists[node->label2_->current_]->delete_node(node);
        delete node;
        return deleted_labels;
    }else
    {
        activeLists[node->label2_->current_]->delete_node(node);
        delete node;
        return 1;
    }
}

int dominance_check(
        std::vector<LabelList*>&    activeLists,
        std::vector<LabelList*>&    oldLists,
        Label*                     label,
        SCIP_Bool                   usedLabels,
        ObjPricerVRP*               pricerData,
        SCIP*                       scip
)
{
    LabelNode* tmp_node;
    LabelNode* next_node;
    int count = 0;

    if(!usedLabels)
    {
        tmp_node = activeLists[label->current_]->head_;
    }else{
        tmp_node = oldLists[label->current_]->head_;
    }

    while (tmp_node != nullptr)
    {
//        tmp_node->label_->print();
        assert(tmp_node != tmp_node->next_);
        next_node = tmp_node->next_;
        if(next_node != nullptr)
            assert(next_node != next_node->next_);
        /* new label gets dominated by an already existing one */
        if(tmp_node->label2_->dominates(label, pricerData, scip))
        {
            assert(count == 0);
            return -1;
        }
        /* new label dominates an already existing one */
        else if (label->dominates(tmp_node->label2_, pricerData, scip))
        {
            /* just delete label if it is non-propagated */
            if(!usedLabels)
            {
                count++;
                activeLists[label->current_]->delete_node(tmp_node);
                delete tmp_node;
            }
                /* delete label and its descendants */
            else
            {
                assert(tmp_node->is_propagated_);
                count += deleteChildren2(activeLists, oldLists, tmp_node, &next_node);
            }
        }
        tmp_node = next_node;
    }

    return count;
}

int dominance_check_enumerate(
        std::vector<LabelList*>&    activeLists,
        std::vector<LabelList*>&    oldLists,
        Label*                     label,
        SCIP_Bool                   usedLabels,
        SCIP*                       scip
)
{
    LabelNode* tmp_node;
    LabelNode* next_node;
    int count = 0;

    if(!usedLabels)
    {
        tmp_node = activeLists[label->current_]->head_;
    }else{
        tmp_node = oldLists[label->current_]->head_;
    }

    while (tmp_node != nullptr)
    {
//        tmp_node->label_->print();
        assert(tmp_node != tmp_node->next_);
        next_node = tmp_node->next_;
        if(next_node != nullptr)
            assert(next_node != next_node->next_);
        /* new label gets dominated by an already existing one */
        if(tmp_node->label2_->dominates_enumerate(label, scip))
        {
            assert(count == 0);
            return -1;
        }
            /* new label dominates an already existing one */
        else if (label->dominates_enumerate(tmp_node->label2_, scip))
        {
            /* just delete label if it is non-propagated */
            if(!usedLabels)
            {
                count++;
                activeLists[label->current_]->delete_node(tmp_node);
                delete tmp_node;
            }
                /* delete label and its descendants */
            else
            {
                assert(tmp_node->is_propagated_);
                count += deleteChildren2(activeLists, oldLists, tmp_node, &next_node);
            }
        }
        tmp_node = next_node;
    }

    return count;
}