"""本地进度；集合以去重后的 JSON 列表保存。"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

DEFAULT_PATH = Path(__file__).resolve().parent.parent / "data" / "progress.json"


def empty_progress():
    return {"chapters": {}, "mistakes": [], "flashcards_mastered": [], "diagnostic": {}}


def load(path=None):
    target = Path(path) if path is not None else DEFAULT_PATH
    if not target.exists():
        return empty_progress()
    try:
        data = json.loads(target.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise ValueError("顶层必须为对象")
        result = empty_progress()
        result.update(data)
        if not isinstance(result["chapters"], dict) or not isinstance(result["diagnostic"], dict):
            raise ValueError("章节与诊断记录必须为对象")
        for key in ("mistakes", "flashcards_mastered"):
            if not isinstance(result[key], list):
                raise ValueError(f"{key} 必须为列表")
        for chapter in result["chapters"].values():
            if not isinstance(chapter, dict):
                raise ValueError("章节记录必须为对象")
            for key in ("code_ex", "quiz_correct_ids", "written"):
                if key in chapter and not isinstance(chapter[key], list):
                    raise ValueError(f"{key} 必须为列表")
        return result
    except (ValueError, TypeError) as exc:
        raise ValueError(f"进度文件格式错误，请先修复或备份 {target}：{exc}") from exc


def save(data, path=None):
    target = Path(path) if path is not None else DEFAULT_PATH
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".progress-", suffix=".json", dir=target.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(data, stream, ensure_ascii=False, indent=2)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def chapter_progress(data, chapter_id):
    return data.get("chapters", {}).get(chapter_id, {
        "theory_done": False, "code_ex": [], "quiz": 0,
        "quiz_correct_ids": [], "written": [],
    })


def _chapter(data, chapter_id):
    return data.setdefault("chapters", {}).setdefault(chapter_id, {
        "theory_done": False, "code_ex": [], "quiz": 0,
        "quiz_correct_ids": [], "written": [],
    })


def _toggle(items, item_id, done):
    items[:] = sorted((set(items) | {item_id}) if done else (set(items) - {item_id}))


def mark_theory(data, chapter_id, done=True):
    _chapter(data, chapter_id)["theory_done"] = bool(done)


def mark_code(data, chapter_id, exercise_id):
    _toggle(_chapter(data, chapter_id).setdefault("code_ex", []), exercise_id, True)


def mark_written(data, chapter_id, exercise_id, done=True):
    _toggle(_chapter(data, chapter_id).setdefault("written", []), exercise_id, done)


def mark_quiz(data, chapter_id, question_id, correct):
    chapter = _chapter(data, chapter_id)
    _toggle(chapter.setdefault("quiz_correct_ids", []), question_id, correct)
    chapter["quiz"] = len(chapter["quiz_correct_ids"])


def set_mistake(data, question_id, included=True):
    _toggle(data.setdefault("mistakes", []), question_id, included)


def set_flashcard(data, card_id, mastered=True):
    _toggle(data.setdefault("flashcards_mastered", []), card_id, mastered)


def set_diagnostic(data, results):
    data["diagnostic"] = dict(results)


def completion(data, chapter):
    """理论、代码、测验、书面题按可完成项目数等权计算。"""
    state = chapter_progress(data, chapter["id"])
    code_ids = {x["id"] for x in chapter["code_exercises"]}
    written_ids = {x["id"] for x in chapter["written"]}
    quiz_ids = {x["id"] for x in chapter["quiz"]}
    total = 1 + len(code_ids) + len(written_ids) + len(quiz_ids)
    done = (int(bool(state.get("theory_done")))
            + len(code_ids.intersection(state.get("code_ex", [])))
            + len(written_ids.intersection(state.get("written", [])))
            + len(quiz_ids.intersection(state.get("quiz_correct_ids", []))))
    return done / total
